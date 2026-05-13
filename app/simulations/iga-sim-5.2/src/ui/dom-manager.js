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
        const muInput = document.getElementById('input-mu');
        const txInput = document.getElementById('input-tx');
        const tolInput = document.getElementById('input-ecsw-tol');
        const tolVal = document.getElementById('ecsw-tol-val');
        if (tolInput) {
            tolInput.oninput = (e) => {
                const exp = parseFloat(e.target.value);
                const val = Math.pow(10, exp);
                tolVal.textContent = val.toExponential(0);
                window.app.params.ecswTol = val;
            };
        }

        const meshToggle = document.getElementById('toggle-ecsw-mesh');
        if (meshToggle) {
            meshToggle.onchange = (e) => {
                window.app.params.showActiveMesh = e.target.checked;
                window.app.requestUpdate();
            };
        }

        rInput.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('r-val').textContent = val.toFixed(2);
            window.app.params.r = val;
            window.app.requestUpdate();
        };

        // Traction Load (Moved to manual trigger)
        txInput.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('tx-val').textContent = val.toFixed(1);
            window.app.params.Tx = val;
            
            // Real-time visual feedback for arrows (no solve)
            if (window.viz && window.app.activePatch) {
                window.viz.renderBCs(window.app.activePatch);
            }
        };

        muInput.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('mu-val').textContent = val.toFixed(2);
            window.app.params.mu = val;
            window.app.requestUpdate();
        };



        const hInput = document.getElementById('input-h');
        const pInput = document.getElementById('input-p');

        hInput.oninput = (e) => {
            document.getElementById('h-val').textContent = e.target.value;
        };

        hInput.onchange = (e) => {
            window.app.params.h = parseInt(e.target.value);
            window.app.requestUpdate();
        };

        pInput.oninput = (e) => {
            document.getElementById('p-val').textContent = parseInt(e.target.value) + 2;
        };

        pInput.onchange = (e) => {
            window.app.params.p = parseInt(e.target.value);
            window.app.requestUpdate();
        };

        // Apply Mesh Button
        document.getElementById('btn-apply-mesh').onclick = () => {
            window.app.requestUpdate();
        };

        // View Mode Toggle
        document.getElementById('view-disp').onclick = () => {
            window.viz.viewMode = 'displacement';
            this.updateViewModeUI();
            window.app.requestUpdate();
        };

        document.getElementById('view-stress').onclick = () => {
            window.viz.viewMode = 'stress';
            this.updateViewModeUI();
            window.app.requestUpdate();
        };

        // Evaluation Mode Toggle
        document.getElementById('mode-fom').onclick = () => {
            window.app.evaluationMode = 'fom';
            this.updateModeUI();
            window.app.requestUpdate();
        };

        document.getElementById('mode-rom').onclick = () => {
            if (window.app.isTrained) {
                window.app.evaluationMode = 'rom';
                this.updateModeUI();
                window.app.requestUpdate();
            }
        };

        document.getElementById('mode-ecsw').onclick = () => {
            if (window.app.isTrained && window.app.ecswData) {
                window.app.evaluationMode = 'ecsw';
                this.updateModeUI();
                window.app.requestUpdate();
            }
        };

        // Train Button
        document.getElementById('btn-train').onclick = () => {
            window.app.runSweep();
        };

        // Toggles
        document.getElementById('toggle-cp').onchange = (e) => {
            if (window.viz) window.viz.setControlNetVisibility(e.target.checked);
        };

        document.getElementById('toggle-bc').onchange = (e) => {
            if (window.viz) window.viz.setBCVisibility(e.target.checked);
        };

        // Tabs
        document.querySelectorAll('.tab').forEach(tab => {
            tab.onclick = () => {
                console.log("[DOM] Tab clicked:", tab.dataset.tab);
                document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
                tab.classList.add('active');
                
                const target = tab.dataset.tab;
                console.log("[DOM] Switching display for:", target);
                document.getElementById('audit-log').style.display = target === 'log' ? 'block' : 'none';
                document.getElementById('reduction-audit').style.display = target === 'reduction' ? 'block' : 'none';

                if (target === 'reduction') {
                    if (window.app && window.app.charts && window.app.charts.energyChart) {
                        window.app.charts.energyChart.resize();
                        window.app.charts.energyChart.update();
                    }
                }
            };
        });
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

    enableROM(enabled) {
        const btn = document.getElementById('mode-rom');
        btn.disabled = !enabled;
        if (enabled) btn.style.color = '#94a3b8';
    }

    enableECSW(enabled) {
        const btn = document.getElementById('mode-ecsw');
        if (btn) {
            btn.disabled = !enabled;
            if (enabled) btn.style.color = '#94a3b8';
        }
    }

    updateModeUI() {
        const mode = window.app.evaluationMode;
        const fomBtn = document.getElementById('mode-fom');
        const romBtn = document.getElementById('mode-rom');
        const ecswBtn = document.getElementById('mode-ecsw');

        [fomBtn, romBtn, ecswBtn].forEach(btn => {
            if (!btn) return;
            btn.style.background = 'transparent';
            btn.style.color = '#64748b';
        });

        const activeBtn = document.getElementById(`mode-${mode}`);
        if (activeBtn) {
            activeBtn.style.background = 'var(--primary)';
            activeBtn.style.color = 'white';
        }
    }

    updateViewModeUI() {
        const isStress = window.viz.viewMode === 'stress';
        const dispBtn = document.getElementById('view-disp');
        const stressBtn = document.getElementById('view-stress');

        if (isStress) {
            stressBtn.style.background = 'var(--primary)';
            stressBtn.style.color = 'white';
            dispBtn.style.background = 'transparent';
            dispBtn.style.color = '#64748b';
        } else {
            dispBtn.style.background = 'var(--primary)';
            dispBtn.style.color = 'white';
            stressBtn.style.background = 'transparent';
            stressBtn.style.color = '#64748b';
        }
    }
}

window.dom = new DOMManager();
