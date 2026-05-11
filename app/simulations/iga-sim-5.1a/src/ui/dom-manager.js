/**
 * Phase 5.1a - DOM Manager
 */

class DOMManager {
    constructor() {
        this.setupListeners();
    }

    setupListeners() {
        // Geometry Sliders
        const rInput = document.getElementById('input-r');
        const l1Input = document.getElementById('input-l1');
        const l2Input = document.getElementById('input-l2');

        rInput.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('r-val').textContent = val.toFixed(2);
            window.app.params.R = val;
            window.app.updateGeometry();
        };

        l1Input.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('l1-val').textContent = val.toFixed(2);
            window.app.params.L1 = val;
            window.app.updateGeometry();
        };

        l2Input.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('l2-val').textContent = val.toFixed(2);
            window.app.params.L2 = val;
            window.app.updateGeometry();
        };

        // Train Button
        document.getElementById('btn-train').onclick = () => {
            window.app.runSweep();
        };
    }

    updateStats(patch, snapCount) {
        const nDofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;
        document.getElementById('stat-dofs').textContent = nDofs;
        document.getElementById('stat-snaps').textContent = snapCount;
    }

    setLoading(msg) {
        document.getElementById('training-status').textContent = msg;
        if (window.app.reporter) window.app.reporter.log('system', msg);
    }

    enableExport(enabled) {
        document.getElementById('btn-export').disabled = !enabled;
    }
}

window.dom = new DOMManager();
