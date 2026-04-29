/**
 * DEIM Online Audit Module
 * Performs a live solve comparison between FOM and DEIM at a target load.
 */
DEIMBenchmarkApp.prototype.runOnlineAudit = async function() {
    console.log("%c--- DEIM ONLINE STRESS TEST ---", "font-weight:bold; color:#ec4899; font-size:1.2em;");
    
    if (!this.isTrained) {
        console.error("Please train DEIM first!");
        return;
    }

    const loadMag = this.loadMag;
    console.log(`Target Load: ${loadMag} N`);
    console.log(`Configuration: k=${this.k}, m=${this.deimM}`);

    // 1. Run FOM (Ground Truth)
    console.log("Step 1: Computing Full Order Solution...");
    const t0 = performance.now();
    const fom = this.solve('fom', loadMag);
    const dt_fom = performance.now() - t0;
    console.log(`✓ FOM Converged in ${fom.result.residualHistory.length} iterations (${dt_fom.toFixed(1)}ms)`);

    // 2. Run DEIM (Reduced Order)
    console.log("Step 2: Computing DEIM Reduced Solution...");
    const t1 = performance.now();
    const deim = this.solve('deim', loadMag);
    const dt_deim = performance.now() - t1;
    console.log(`✓ DEIM Converged in ${deim.result.residualHistory.length} iterations (${dt_deim.toFixed(1)}ms)`);

    // 3. Error Analysis
    let num = 0, den = 0;
    for (let i = 0; i < fom.result.u.length; i++) {
        num += (fom.result.u[i] - deim.result.u[i])**2;
        den += fom.result.u[i]**2;
    }
    const relError = Math.sqrt(num/den);
    const speedup = dt_fom / dt_deim;

    console.log("%c--- RESULTS ---", "font-weight:bold;");
    console.log(`L2 Relative Error : ${(relError * 100).toFixed(6)}%`);
    console.log(`Runtime Speedup   : ${speedup.toFixed(2)}x`);

    if (relError > 0.05) {
        console.warn("⚠ WARNING: High Online Error detected (>5%).");
        console.log("Potential Causes:");
        console.log(" - Load Mag is outside training range.");
        console.log(" - k modes (basis) is too small for this deformation level.");
        console.log(" - Interpolation points m is insufficient for this load.");
    } else {
        console.log("%c✓ DEIM Online Accuracy Verified", "color:green; font-weight:bold;");
    }
};
