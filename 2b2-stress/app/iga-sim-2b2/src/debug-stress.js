// Enhanced stress diagnostic - checks forces, displacement consistency, and stress at key points
function debugStress2() {
    console.log("=== ENHANCED STRESS DIAGNOSTIC ===");
    console.log("E:", targetState.E, "nu:", targetState.nu, "Tx:", targetState.load);
    
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    const nDofs = nU * nV * 2;
    
    // 1. Check applied forces
    const forces = solver.calculateNodalTraction(patch, targetState.load, 'right');
    console.log("\n=== APPLIED FORCES (from Kirsch traction) ===");
    let totalFx = 0, totalFy = 0;
    for (let i = 0; i < nU; i++) {
        for (let j = 0; j < nV; j++) {
            const idx = (i * nV + j) * 2;
            const fx = forces[idx];
            const fy = forces[idx + 1];
            if (Math.abs(fx) > 0.001 || Math.abs(fy) > 0.001) {
                const cp = patch.controlPoints[i][j];
                console.log(`  CP[${i}][${j}] at (${cp.x.toFixed(2)}, ${cp.y.toFixed(2)}): fx=${fx.toFixed(3)}, fy=${fy.toFixed(3)}`);
                totalFx += fx;
                totalFy += fy;
            }
        }
    }
    console.log(`  TOTAL: Fx=${totalFx.toFixed(3)}, Fy=${totalFy.toFixed(3)}`);
    
    // 2. Check physical displacement at sample points (not CP displacement)
    console.log("\n=== PHYSICAL DISPLACEMENT (interpolated) ===");
    const samplePts = [
        {u: 0.0, v: 0.0}, {u: 0.0, v: 0.5}, {u: 0.0, v: 1.0},
        {u: 0.5, v: 0.0}, {u: 0.5, v: 0.5},
        {u: 1.0, v: 0.0}, {u: 1.0, v: 0.5}, {u: 1.0, v: 1.0}
    ];
    for (const sp of samplePts) {
        const state = engine.getSurfaceState(patch, sp.u, sp.v);
        const interp = interpolateDisplacement(sp.u, sp.v, state.denominator);
        console.log(`  (u=${sp.u}, v=${sp.v}) -> phys (${state.position.x.toFixed(3)}, ${state.position.y.toFixed(3)}) -> disp (${interp.x.toFixed(6)}, ${interp.y.toFixed(6)})`);
    }
    
    // 3. Verify stress at mid-field point where exact solution is simple
    // At (L, 0), the exact solution should be close to far-field: σ_xx ≈ Tx
    console.log("\n=== STRESS VERIFICATION ===");
    const verifyPts = [
        {u: 0.0, v: 0.5, label: "mid-bottom (≈2.5, 0)"},
        {u: 0.0, v: 0.9, label: "near-bottom-outer (≈3.7, 0)"},
        {u: 1.0, v: 0.5, label: "mid-left (≈0, 2.5)"},
    ];
    for (const vp of verifyPts) {
        const pos = engine.evaluateSurface(patch, vp.u, vp.v);
        const sN = solver.getNumericalStress(patch, analysisData.u, vp.u, vp.v, targetState.E, targetState.nu);
        const sE = solver.getExactStress(pos.x, pos.y, 1.0, targetState.load);
        const ratio_xx = sE.sxx !== 0 ? (sN.sxx / sE.sxx * 100).toFixed(1) : 'N/A';
        console.log(`  ${vp.label}: Num_sxx=${sN.sxx.toFixed(2)}, Exact_sxx=${sE.sxx.toFixed(2)}, ratio=${ratio_xx}%`);
    }
    
    // 4. Check stiffness matrix condition
    console.log("\n=== STIFFNESS MATRIX DIAGNOSTIC ===");
    const K = solver.assembleStiffness(patch);
    let maxK = 0, minNonZeroK = Infinity;
    for (let i = 0; i < nDofs; i++) {
        for (let j = 0; j < nDofs; j++) {
            const val = Math.abs(K[i][j]);
            if (val > maxK) maxK = val;
            if (val > 1e-15 && val < minNonZeroK) minNonZeroK = val;
        }
    }
    console.log(`  Max |K_ij| = ${maxK.toExponential(4)}`);
    console.log(`  Min nonzero |K_ij| = ${minNonZeroK.toExponential(4)}`);
    console.log(`  Condition ratio ≈ ${(maxK/minNonZeroK).toExponential(4)}`);
    console.log(`  Diagonal K[0][0] = ${K[0][0].toFixed(4)}, K[1][1] = ${K[1][1].toFixed(4)}`);
}

window.debugStress2 = debugStress2;
console.log("Enhanced debug loaded. Run debugStress2() after Analysis in Stress mode.");
