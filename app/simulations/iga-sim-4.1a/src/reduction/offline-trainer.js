/**
 * Phase 4.1 Offline Trainer
 * Handles POD extraction, Reduced Matrix generation, and JSON Export.
 */

class OfflineTrainer {
    constructor(app) {
        this.app = app;
        this.package = {
            metadata: { date: new Date().toISOString(), version: "4.1b-online" },
            phi: null,
            ecsw: { weights: [], indices: [] },
            reducedMatrices: { M_r: null, C_r: null },
            mesh: { h: 0, p: 1 },
            k: 0,
            snapshots: 0,
            pod_energy: 0
        };
    }

    /**
     * Compute POD Basis from displacement snapshots
     */
    computePOD(snapshots, k) {
        console.log(`[Offline] Computing POD on ${snapshots.length} snapshots...`);
        const n = snapshots[0].length;
        const S = new mlMatrix.Matrix(n, snapshots.length);
        
        for (let j = 0; j < snapshots.length; j++) {
            for (let i = 0; i < n; i++) {
                S.set(i, j, snapshots[j][i]);
            }
        }

        // SVD: S = U * Sigma * V^T
        const svd = new mlMatrix.SingularValueDecomposition(S, { autoTranspose: true });
        const U = svd.leftSingularVectors;
        const s = svd.diagonal;

        // Energy check
        const totalVar = s.reduce((a, b) => a + b*b, 0);
        const energyTrace = [];
        
        if (totalVar < 1e-20) {
            console.warn("[Offline] Near-zero variance detected. Defaulting energy to 100%.");
            for (let i = 0; i < s.length; i++) energyTrace.push(1.0);
        } else {
            let cumVar = 0;
            for (let i = 0; i < Math.min(s.length, 30); i++) {
                cumVar += s[i]*s[i];
                energyTrace.push(cumVar / totalVar);
            }
        }

        // Truncate to k modes (ensure k doesn't exceed available columns)
        const actualK = Math.min(k, U.columns, s.length);
        const Phi = U.subMatrix(0, n - 1, 0, actualK - 1);
        this.package.phi = Phi.to2DArray();
        this.package.pod_energy = energyTrace[actualK - 1];
        this.package.k = actualK;
        this.package.snapshots = snapshots.length;
        
        console.log(`[Offline] POD Complete. Extracted ${actualK} modes. Energy: ${(this.package.pod_energy * 100).toFixed(4)}%`);
        return { Phi, energyTrace };
    }

    /**
     * Phase 3: ECSW Sparse Sampling
     */
    async trainECSW(snapshots, Phi, tol = 1e-4) {
        console.log(`[Offline] Starting ECSW Training with tol=${tol}...`);
        
        // 1. Subsample snapshots equidistantly to speed up training
        // If we have 100 snapshots, training on all might be slow.
        // We pick a subset that covers the entire force range equidistantly.
        const maxTrainSnaps = 20;
        const trainSnapshots = [];
        if (snapshots.length <= maxTrainSnaps) {
            snapshots.forEach(s => trainSnapshots.push(s));
        } else {
            const step = (snapshots.length - 1) / (maxTrainSnaps - 1);
            for (let i = 0; i < maxTrainSnaps; i++) {
                trainSnapshots.push(snapshots[Math.round(i * step)]);
            }
        }
        console.log(`[Offline] Subsampled ${trainSnapshots.length}/${snapshots.length} snapshots for ECSW.`);

        // 2. Initialize and Run ECSW
        const ecsw = new window.ECSWEngine();
        const romDummy = { Phi: Phi }; 
        
        const result = await ecsw.train(this.app.fom, romDummy, this.app.patch, trainSnapshots, tol);
        
        // 3. Store in package
        this.package.ecsw = {
            indices: ecsw.sampleElements.map(el => el.i * ( [...new Set(this.app.patch.V)].length - 1) + el.j),
            weights: ecsw.sampleElements.map(el => el.weight)
        };
        
        // Store mesh info
        this.package.mesh = {
            h: parseInt(document.getElementById('input-h').value),
            p: parseInt(document.getElementById('input-p').value)
        };

        return { 
            indices: this.package.ecsw.indices, 
            weights: this.package.ecsw.weights 
        };
    }

    /**
     * Phase 4: Project Global Matrices
     */
    computeReducedMatrices(Phi) {
        const M = new mlMatrix.Matrix(this.app.dyn.M);
        const PhiT = Phi.transpose();
        const Mr = PhiT.mmul(M).mmul(Phi);
        const Cr = Mr.clone().mul(this.app.dyn.alpha); 
        
        this.package.reducedMatrices.M_r = Mr.to2DArray();
        this.package.reducedMatrices.C_r = Cr.to2DArray();
        this.package.dofs = this.app.dyn.dofs;
    }

    exportPackage() {
        const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(this.package, null, 2));
        const downloadAnchorNode = document.createElement('a');
        downloadAnchorNode.setAttribute("href", dataStr);
        downloadAnchorNode.setAttribute("download", "rom_package_4.1b.json");
        document.body.appendChild(downloadAnchorNode);
        downloadAnchorNode.click();
        downloadAnchorNode.remove();
    }

    verifyPackage() {
        const pkg = this.package;
        const report = {
            isValid: true,
            checks: []
        };

        // 1. Basis Existence & Dimensions
        if (pkg.phi) {
            const Phi = new mlMatrix.Matrix(pkg.phi);
            report.checks.push({ name: "Basis Dimensions", status: "OK", detail: `${Phi.rows} DOFs x ${Phi.columns} Modes` });
            
            // 2. Orthogonality Check: Phi^T * Phi should be I
            const PhiT = Phi.transpose();
            const I_check = PhiT.mmul(Phi);
            let orthError = 0;
            for(let i=0; i<I_check.rows; i++) {
                for(let j=0; j<I_check.columns; j++) {
                    const expected = (i === j) ? 1.0 : 0.0;
                    orthError = Math.max(orthError, Math.abs(I_check.get(i,j) - expected));
                }
            }
            const orthStatus = orthError < 1e-10 ? "OK" : "Warning";
            report.checks.push({ name: "Basis Orthogonality", status: orthStatus, detail: `Max Error: ${orthError.toExponential(2)}` });
        } else {
            report.isValid = false;
            report.checks.push({ name: "Basis Existence", status: "FAIL", detail: "Missing POD Basis" });
        }

        // 3. Energy Check
        const energy = (pkg.pod_energy || 0) * 100;
        const energyStatus = energy > 99.9 ? "OK" : "Warning";
        report.checks.push({ name: "Energy Capture", status: energyStatus, detail: `${energy.toFixed(4)}%` });

        // 4. Reduced Matrices
        if (pkg.reducedMatrices.M_r) {
            const Mr = new mlMatrix.Matrix(pkg.reducedMatrices.M_r);
            report.checks.push({ name: "Reduced Mass Matrix", status: "OK", detail: `${Mr.rows}x${Mr.columns} (Dense)` });
        } else {
            report.isValid = false;
            report.checks.push({ name: "Matrix Integrity", status: "FAIL", detail: "Missing Reduced Matrices" });
        }

        return report;
    }
}
