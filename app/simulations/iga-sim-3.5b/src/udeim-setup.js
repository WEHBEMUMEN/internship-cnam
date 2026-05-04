/**
 * Phase 3.4b: U-DEIM Setup/Projection
 */

UDEIMEngine.prototype.precomputeReducedTangent = function(fomSolver, romEngine, patch, snapU, bcs = []) {
    const { Matrix } = window.mlMatrix;
    const Phi = romEngine.Phi;
    const k = Phi.columns;
    const PhiT = Phi.transpose();

    const Kt_red_avg = Array.from({ length: k }, () => new Float64Array(k));
    const nSnaps = snapU.length;

    for (let s = 0; s < nSnaps; s++) {
        const Kt_full = fomSolver.calculateTangentStiffness(patch, snapU[s]);
        fomSolver.applyPenaltyConstraints(Kt_full, null, snapU[s], patch, bcs);
        const Kt_mat = new Matrix(Kt_full.map(r => Array.from(r)));
        const Kt_red = PhiT.mmul(Kt_mat).mmul(Phi).to2DArray();
        for (let i = 0; i < k; i++)
            for (let j = 0; j < k; j++)
                Kt_red_avg[i][j] += Kt_red[i][j] / nSnaps;
    }
    this.Kt_red_ref = Kt_red_avg;
};
