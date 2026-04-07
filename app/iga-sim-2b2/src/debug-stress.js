// Debug stress recovery - run this in the browser console after solving
// This will diagnose why the L2 error is 93%

function debugStress() {
    console.log("=== STRESS DIAGNOSTIC ===");
    console.log("Patch:", patch.p, patch.q, "nU:", patch.controlPoints.length, "nV:", patch.controlPoints[0].length);
    console.log("E:", targetState.E, "nu:", targetState.nu, "Tx:", targetState.load);
    
    // Test point: mid-hole boundary at theta=90 (top of hole)
    // In parametric space: u=1 (left edge, 90 deg), v=0 (hole boundary)
    const testPoints = [
        { u: 0.0, v: 0.0, label: "u=0,v=0 (R,0) bottom-hole" },
        { u: 1.0, v: 0.0, label: "u=1,v=0 (0,R) left-hole" },
        { u: 0.5, v: 0.0, label: "u=0.5,v=0 (hole at 45deg)" },
        { u: 0.0, v: 0.5, label: "u=0,v=0.5 (bottom mid-radial)" },
        { u: 0.5, v: 0.5, label: "u=0.5,v=0.5 (center)" },
        { u: 0.0, v: 1.0, label: "u=0,v=1 (L,0) bottom-outer" },
    ];
    
    for (const tp of testPoints) {
        const pos = engine.evaluateSurface(patch, tp.u, tp.v);
        const sN = solver.getNumericalStress(patch, analysisData.u, tp.u, tp.v, targetState.E, targetState.nu);
        const sE = solver.getExactStress(pos.x, pos.y, 1.0, targetState.load);
        
        console.log(`\n--- ${tp.label} ---`);
        console.log(`  Physical: (${pos.x.toFixed(3)}, ${pos.y.toFixed(3)})`);
        console.log(`  Numerical: sxx=${sN.sxx.toFixed(2)}, syy=${sN.syy.toFixed(2)}, sxy=${sN.sxy.toFixed(2)}, vM=${sN.vonMises.toFixed(2)}`);
        console.log(`  Exact:     sxx=${sE.sxx.toFixed(2)}, syy=${sE.syy.toFixed(2)}, sxy=${sE.sxy.toFixed(2)}`);
        console.log(`  Error:     dxx=${(sN.sxx - sE.sxx).toFixed(2)}, dyy=${(sN.syy - sE.syy).toFixed(2)}, dxy=${(sN.sxy - sE.sxy).toFixed(2)}`);
    }
    
    // Also check the displacement at a few known points
    console.log("\n=== DISPLACEMENT CHECK ===");
    const nV = patch.controlPoints[0].length;
    for (let i = 0; i < patch.controlPoints.length; i++) {
        for (let j = 0; j < nV; j++) {
            const idx = (i * nV + j) * 2;
            const cp = patch.controlPoints[i][j];
            if (Math.abs(analysisData.u[idx]) > 0.001 || Math.abs(analysisData.u[idx+1]) > 0.001) {
                console.log(`  CP[${i}][${j}] at (${cp.x.toFixed(2)}, ${cp.y.toFixed(2)}): ux=${analysisData.u[idx].toFixed(5)}, uy=${analysisData.u[idx+1].toFixed(5)}`);
            }
        }
    }
}

window.debugStress = debugStress;
console.log("Debug loaded. Run debugStress() after switching to Stress mode and pressing Analysis.");
