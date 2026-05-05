/**
 * Phase 3.4b: U-DEIM Trainer
 */

UDEIMEngine.prototype.train = async function(fomSolver, romEngine, patch, snapshotDisplacements, m, kf = m, bcs = []) {
    this.m = m;
    this.kf = Math.min(kf, m);
    this.bcs = bcs; // Store for audit/offline use

    
    const nSnaps = snapshotDisplacements.length;
    const { U, V, controlPoints, p, q } = patch;
    const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
    
    const elements = [];
    const nV = controlPoints[0].length;
    for (let i = 0; i < uniqueU.length - 1; i++) {
        if (uniqueU[i+1] - uniqueU[i] < 1e-10) continue;
        for (let j = 0; j < uniqueV.length - 1; j++) {
            if (uniqueV[j+1] - uniqueV[j] < 1e-10) continue;
            elements.push({ 
                index: elements.length,
                i, j, 
                uMin: uniqueU[i], uMax: uniqueU[i+1],
                vMin: uniqueV[j], vMax: uniqueV[j+1]
            });
        }
    }
    this.numElements = elements.length;
    this.numLocalDofs = (p + 1) * (q + 1) * 2;
    const N_u = this.numElements * this.numLocalDofs; // Size of unassembled force vector
    
    // 1. Build map from (element, localDOF) -> globalDOF
    this.elementDofMap = [];
    const nU = patch.controlPoints.length;
    
    for (let e = 0; e < this.numElements; e++) {
        const map = new Int32Array(this.numLocalDofs);
        const el = elements[e];
        const uMid = (el.uMin + el.uMax) / 2.0;
        const vMid = (el.vMin + el.vMax) / 2.0;
        
        const kU = fomSolver.engine.findSpan(nU - 1, p, uMid, U);
        const kV = fomSolver.engine.findSpan(nV - 1, q, vMid, V);
        
        let ld = 0;
        for (let ii = kU - p; ii <= kU; ii++) {
            for (let jj = kV - q; jj <= kV; jj++) {
                const cpIdx = ii * nV + jj;
                map[ld++] = cpIdx * 2;
                map[ld++] = cpIdx * 2 + 1;
            }
        }
        this.elementDofMap.push(map);
    }

    // 2. Compute Unassembled Force Snapshots
    const S_u = [];
    const gRule = window.GaussQuadrature2D.getPoints(Math.max(p, q) + 1);
    
    console.log(`U-DEIM: Computing unassembled forces for ${nSnaps} snapshots...`);
    for (let s = 0; s < nSnaps; s++) {
        const f_u = new Float64Array(N_u);
        const u_disp = snapshotDisplacements[s];
        
        let normSq = 0;
        for (let e = 0; e < this.numElements; e++) {
            const el = elements[e];
            const f_e = this._computeLocalForce(fomSolver, patch, el, u_disp, gRule);
            for (let ld = 0; ld < this.numLocalDofs; ld++) {
                const val = f_e[ld];
                f_u[e * this.numLocalDofs + ld] = val;
                normSq += val * val;
            }
        }
        
        // Normalize snapshots to avoid high-load dominance in SVD
        const norm = Math.sqrt(normSq);
        if (norm > 1e-12) {
            for (let i = 0; i < N_u; i++) f_u[i] /= norm;
        }
        S_u.push(f_u);
    }

    // 3. POD on Unassembled Forces
    console.log(`U-DEIM: Performing SVD on ${N_u}x${nSnaps} snapshot matrix...`);
    const { Matrix, SVD } = window.mlMatrix;
    const S_mat = new Matrix(S_u.map(s => Array.from(s))).transpose();
    const svd = new SVD(S_mat, { computeLeftSingularVectors: true, computeRightSingularVectors: false, autoTranspose: true });
    const U_greedy = svd.leftSingularVectors;
    
    this.m = Math.min(this.m, U_greedy.columns);
    this.kf = Math.min(this.kf, this.m);
    this.U_f = U_greedy.subMatrix(0, N_u - 1, 0, this.kf - 1);
    
    // 4. Q-DEIM: QR Factorization with Column Pivoting on U_f^T
    // This corresponds to performing pivoted Gram-Schmidt directly on the rows of U_f.
    const A = Array.from({ length: N_u }, (_, i) => {
        const row = new Float64Array(this.m);
        for (let j = 0; j < this.m; j++) row[j] = U_greedy.get(i, j);
        return row;
    });
    
    // Keep track of the squared 2-norm of each row (unassembled DOF)
    const norms = new Float64Array(N_u);
    for (let i = 0; i < N_u; i++) {
        let sum = 0;
        for (let j = 0; j < this.m; j++) sum += A[i][j] * A[i][j];
        norms[i] = sum;
    }
    
    const constrainedDofs_init = new Set();
    const nV_global_init = patch.controlPoints[0].length;
    if (Array.isArray(this.bcs)) {
        this.bcs.forEach(bc => {
            const base = (bc.i * nV_global_init + bc.j) * 2;
            if (bc.axis === 'x' || bc.axis === 'both') constrainedDofs_init.add(base);
            if (bc.axis === 'y' || bc.axis === 'both') constrainedDofs_init.add(base + 1);
        });
    }

    const indices = [];
    this.history = [];

    for (let l = 0; l < this.m; l++) {
        // Find row with max residual norm (Pivoting)
        let maxNorm = -1, pivot = -1;
        for (let i = 0; i < N_u; i++) {
            // Skip constrained DOFs
            if (constrainedDofs_init.has(this.elementDofMap[Math.floor(i / this.numLocalDofs)][i % this.numLocalDofs])) continue;
            
            if (norms[i] > maxNorm && !indices.includes(i)) {
                maxNorm = norms[i];
                pivot = i;
            }
        }

        if (pivot === -1 || maxNorm < 1e-22) {
            console.warn(`U-DEIM: Residual dropped to noise floor at step ${l+1}. Truncating.`);
            this.m = l;
            this.kf = Math.min(this.kf, this.m);
            break;
        }

        indices.push(pivot);
        this.history.push({
            step: l + 1,
            point: pivot,
            maxVal: Math.sqrt(maxNorm),
            residual: Array.from(norms).map(n => Math.sqrt(Math.max(0, n)))
        });

        // Orthogonalize remaining rows against the pivot row
        const pivotRow = A[pivot];
        const pNorm = Math.sqrt(maxNorm);
        for (let j = 0; j < this.m; j++) pivotRow[j] /= pNorm; 

        for (let i = 0; i < N_u; i++) {
            if (indices.includes(i)) continue;
            
            // Dot product
            let dot = 0;
            for (let j = 0; j < this.m; j++) dot += A[i][j] * pivotRow[j];
            
            // Subtract projection and recompute norm
            let newNorm = 0;
            for (let j = 0; j < this.m; j++) {
                A[i][j] -= dot * pivotRow[j];
                newNorm += A[i][j] * A[i][j];
            }
            norms[i] = newNorm;
        }
    }
    this.indices = indices;

    // Map indices back to specific elements and local DOFs
    this.sampledDofs = indices.map(idx => ({
        e: Math.floor(idx / this.numLocalDofs),
        ld: idx % this.numLocalDofs
    }));
    
    // Find unique active elements
    const uniqueElIndices = [...new Set(this.sampledDofs.map(s => s.e))];
    this.activeElements = uniqueElIndices.map(eIdx => elements[eIdx]);

    // 5. Precompute Projection Matrix M
    // M = Phi^T * A * U_f * (P^T U_f)^+
    console.log("U-DEIM: Precomputing projection operators...");
    
    // P^T U_f (Interpolation matrix)
    const PtU_arr = Array.from({ length: this.m }, () => new Float64Array(this.kf));
    for (let i = 0; i < this.m; i++) {
        for (let j = 0; j < this.kf; j++) PtU_arr[i][j] = this.U_f.get(this.indices[i], j);
    }
    const PtU_mat = new Matrix(PtU_arr.map(r => Array.from(r)));
    
    // Robust Pseudo-Inverse using SVD ( (P^T U_f)^+ )
    const svd_interp = new SVD(PtU_mat, { computeLeftSingularVectors: true, computeRightSingularVectors: true });
    const U_interp = svd_interp.leftSingularVectors;
    const S_interp = svd_interp.diagonal;
    const V_interp = svd_interp.rightSingularVectors;
    
    // S_inv = diag(1/sigma) with threshold
    const S_inv_mat = Matrix.zeros(V_interp.columns, U_interp.columns);
    const threshold = Math.max(...S_interp) * 1e-12;
    for (let i = 0; i < S_interp.length; i++) {
        if (S_interp[i] > threshold) S_inv_mat.set(i, i, 1.0 / S_interp[i]);
    }
    const PtU_pinv = V_interp.mmul(S_inv_mat).mmul(U_interp.transpose());
    this.PtU_pinv = PtU_pinv.to2DArray();

    // Audit conditioning
    const cond = S_interp[0] / Math.max(S_interp[S_interp.length-1], 1e-30);
    console.log(`U-DEIM: Interpolation Matrix Condition Number: ${cond.toExponential(2)} (${cond < 1e6 ? 'HEALTHY' : 'DANGEROUS'})`);


    // A * U_f (Assemble U_f into global space)
    const Phi = romEngine.Phi;
    const nGlobalDofs = Phi.rows;
    
    const A_Uf = new Matrix(Array.from({length: nGlobalDofs}, () => new Array(this.kf).fill(0)));
    for (let j = 0; j < this.kf; j++) {
        for (let e = 0; e < this.numElements; e++) {
            const map = this.elementDofMap[e];
            for (let ld = 0; ld < this.numLocalDofs; ld++) {
                const globalDof = map[ld];
                const val = this.U_f.get(e * this.numLocalDofs + ld, j);
                A_Uf.set(globalDof, j, A_Uf.get(globalDof, j) + val);
            }
        }
    }
    
    // Phi^T * (A * U_f) * PtU_pinv
    const PhiT_A_Uf = Phi.transpose().mmul(A_Uf);
    const M_matrix = PhiT_A_Uf.mmul(new Matrix(this.PtU_pinv));
    this.M = M_matrix.to2DArray();
    
    // Pre-compute Reduced Tangent
    this.precomputeReducedTangent(fomSolver, romEngine, patch, snapshotDisplacements, bcs);



    console.log(`U-DEIM Trained: ${uniqueElIndices.length} elements selected.`);
    return {
        m: this.m,
        history: this.history,
        elementCount: uniqueElIndices.length,
        totalElements: this.numElements
    };
};

UDEIMEngine.prototype._computeLocalForce = function(fomSolver, patch, el, u_disp, gRule) {
    const { uMin, uMax, vMin, vMax } = el;
    const f_e = new Float64Array(this.numLocalDofs);
    let ld = 0;

    for (let gu = 0; gu < gRule.points.length; gu++) {
        const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
        const wu = gRule.weights[gu] * (uMax - uMin) / 2;
        for (let gv = 0; gv < gRule.points.length; gv++) {
            const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
            const wv = gRule.weights[gv] * (vMax - vMin) / 2;

            const deriv = fomSolver.engine.getSurfaceDerivatives(patch, u, v);
            const { grads: B_param, detJ, activeIndices } = fomSolver.getBParametric(patch, u, v, deriv);
            
            let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
            for (let a = 0; a < activeIndices.length; a++) {
                const k = activeIndices[a];
                dudx += B_param[k][0] * u_disp[k * 2];
                dudy += B_param[k][1] * u_disp[k * 2];
                dvdx += B_param[k][0] * u_disp[k * 2 + 1];
                dvdy += B_param[k][1] * u_disp[k * 2 + 1];
            }

            const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
            const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
            const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
            
            const D = fomSolver.getPlaneStressD();
            const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
            const Syy = D[1][0]*Exx + D[1][1]*Eyy;
            const Sxy = D[2][2]*Exy2;

            const factor = detJ * wu * wv * fomSolver.thickness;

            for (let a = 0; a < activeIndices.length; a++) {
                const k = activeIndices[a];
                const dRdx = B_param[k][0], dRdy = B_param[k][1];
                
                const bexx_u = (1 + dudx) * dRdx;
                const bexx_v = (dvdx) * dRdx;
                const beyy_u = (dudy) * dRdy;
                const beyy_v = (1 + dvdy) * dRdy;
                const bexy_u = (1 + dudx)*dRdy + dudy*dRdx;
                const bexy_v = (1 + dvdy)*dRdx + dvdx*dRdy;

                f_e[a * 2] += (bexx_u * Sxx + beyy_u * Syy + bexy_u * Sxy) * factor;
                f_e[a * 2 + 1] += (bexx_v * Sxx + beyy_v * Syy + bexy_v * Sxy) * factor;
            }
        }
    }
    return f_e;
};
