/**
 * DEIM Benchmark — UI Interactions
 */

DEIMBenchmarkApp.prototype.initUI = function() {
    document.getElementById('btn-train').onclick = () => this.trainAll();
    document.getElementById('btn-compare').onclick = () => this.runComparison();
    document.getElementById('btn-audit').onclick = () => this.auditMath();
    document.getElementById('btn-explorer').onclick = () => {
        const explorerTab = document.querySelector('.chart-tab[data-chart="explorer"]');
        if (explorerTab) explorerTab.click();
    };

    // Phase 3.5 Specific Sliders
    const samplesInput = document.getElementById('input-samples');
    const maxForceInput = document.getElementById('input-max-force');
    if (samplesInput) {
        samplesInput.oninput = (e) => {
            document.getElementById('samples-val').textContent = e.target.value;
        };
    }
    if (maxForceInput) {
        maxForceInput.oninput = (e) => {
            document.getElementById('max-force-val').textContent = e.target.value;
            const physicsLoadInput = document.getElementById('input-load');
            if (physicsLoadInput) physicsLoadInput.max = e.target.value;
        };
    }
    
    document.getElementById('input-mesh').onchange = (e) => {
        this.meshLevel = parseInt(e.target.value);
        document.getElementById('mesh-val').textContent = this.meshLevel;
        this.loadBenchmark();
    };

    document.getElementById('input-load').oninput = (e) => {
        this.loadMag = parseFloat(e.target.value);
        document.getElementById('load-val').textContent = this.loadMag;
        this.updatePhysics();
    };

    document.getElementById('input-k').oninput = (e) => {
        this.k = parseInt(e.target.value);
        document.getElementById('k-val').textContent = this.k;
    };

    document.getElementById('input-m').oninput = (e) => {
        this.deimM = parseInt(e.target.value);
        document.getElementById('m-val').textContent = this.deimM;
    };

    document.querySelectorAll('[data-method]').forEach(btn => {
        btn.onclick = () => {
            document.querySelectorAll('[data-method]').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            this.method = btn.dataset.method;
            this.updatePhysics();
        };
    });

    document.getElementById('input-show-dofs').onchange = (e) => {
        this.showDOFs = e.target.checked;
        this.updateMesh(this.lastResult ? this.lastResult.u : null);
        this._render();
    };

    document.getElementById('input-load-type').onchange = (e) => {
        this.loadType = e.target.value;
        this.updatePhysics();
    };

    // Explorer Nav
    document.getElementById('btn-exp-prev').onclick = () => {
        if (this.explorerStep > 1) {
            this.explorerStep--;
            this.runDEIMExplorer(this.explorerStep);
        }
    };
    document.getElementById('btn-exp-next').onclick = () => {
        if (this.explorerStep < this.deimEngine.history.length) {
            this.explorerStep++;
            this.runDEIMExplorer(this.explorerStep);
        }
    };
};

// --- Entry Point ---
window.onload = () => {
    window.app = new DEIMBenchmarkApp();
};
