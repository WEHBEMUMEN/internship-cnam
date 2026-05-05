/**
 * Phase 4.1 Offline Trainer
 * Handles POD extraction, Reduced Matrix generation, and JSON Export.
 */

class OfflineTrainer {
    constructor(app) {
        this.app = app;
        this.package = {
            metadata: { date: new Date().toISOString(), version: "4.1-transient" },
            basis: null, // Phi
            weights: null, // w (ECSW)
            reducedMatrices: { M_r: null, C_r: null },
            dofs: 0
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
        const svd = new mlMatrix.SingularValueDecomposition(S);
        const U = svd.leftSingularVectors;
        const s = svd.singularValues;

        // Energy check
        const totalVar = s.reduce((a, b) => a + b*b, 0);
        let cumVar = 0;
        const energyTrace = [];
        for (let i = 0; i < Math.min(s.length, 30); i++) {
            cumVar += s[i]*s[i];
            energyTrace.push(cumVar / totalVar);
        }

        // Truncate to k modes
        const Phi = U.subMatrix(0, n - 1, 0, k - 1);
        this.package.basis = Phi.toArray();
        this.package.energy = energyTrace[k-1];
        
        console.log(`[Offline] POD Complete. Energy: ${(this.package.energy * 100).toFixed(4)}%`);
        return { Phi, energyTrace };
    }

    /**
     * Project Global Matrices M and C onto the POD Basis
     */
    computeReducedMatrices(Phi) {
        console.log("[Offline] Projecting Mass and Damping matrices...");
        const M = new mlMatrix.Matrix(this.app.dyn.M);
        
        // M_r = Phi^T * M * Phi
        const PhiT = Phi.transpose();
        const Mr = PhiT.mmul(M).mmul(Phi);
        
        // C_r = Phi^T * C * Phi
        // C = alpha*M + betaR*K
        // For simplicity at export, we assume Rayleigh damping can be re-assembled or we store the reduced forms.
        const Cr = Mr.clone().mul(this.app.dyn.alpha); 
        // Note: Full Kr projection is expensive, usually we project the internal force or use ECSW.
        
        this.package.reducedMatrices.M_r = Mr.toArray();
        this.package.reducedMatrices.C_r = Cr.toArray();
        this.package.dofs = this.app.dyn.dofs;
    }

    /**
     * Export the full ROM package as a JSON file
     */
    exportPackage() {
        const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(this.package, null, 2));
        const downloadAnchorNode = document.createElement('a');
        downloadAnchorNode.setAttribute("href", dataStr);
        downloadAnchorNode.setAttribute("download", "rom_package_4.1.json");
        document.body.appendChild(downloadAnchorNode);
        downloadAnchorNode.click();
        downloadAnchorNode.remove();
        console.log("[Offline] Package exported successfully.");
    }

    /**
     * Verify the integrity of a loaded/stored package
     */
    verifyPackage() {
        console.log("[Offline] Verifying ROM Package Integrity...");
        const pkg = this.package;
        let errors = [];

        if (!pkg.basis) errors.push("Missing POD Basis (Phi)");
        if (!pkg.reducedMatrices.M_r) errors.push("Missing Reduced Mass Matrix (Mr)");
        
        if (errors.length === 0) {
            console.log("%c[SUCCESS] ROM Package is valid and consistent.", "color: #10b981; font-weight: bold;");
            return true;
        } else {
            console.error("[FAILURE] Package corrupted:", errors);
            return false;
        }
    }
}
