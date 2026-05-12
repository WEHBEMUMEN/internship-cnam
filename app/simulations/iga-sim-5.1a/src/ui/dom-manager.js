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
        const txInput = document.getElementById('input-tx');

        rInput.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('r-val').textContent = val.toFixed(2);
            window.app.params.R = val;
            window.app.updateGeometry();
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

        const hInput = document.getElementById('input-h');
        const pInput = document.getElementById('input-p');

        hInput.oninput = (e) => {
            document.getElementById('h-val').textContent = e.target.value;
        };

        hInput.onchange = (e) => {
            window.app.params.h = parseInt(e.target.value);
            // Removed automatic updateGeometry
        };

        pInput.oninput = (e) => {
            document.getElementById('p-val').textContent = parseInt(e.target.value) + 2;
        };

        pInput.onchange = (e) => {
            window.app.params.p = parseInt(e.target.value);
            // Removed automatic updateGeometry
        };

        // Apply Mesh Button
        document.getElementById('btn-apply-mesh').onclick = () => {
            window.app.updateGeometry();
        };

        // View Mode Toggle
        document.getElementById('view-disp').onclick = () => {
            window.viz.viewMode = 'displacement';
            this.updateViewModeUI();
            window.app.updateGeometry();
        };

        document.getElementById('view-stress').onclick = () => {
            window.viz.viewMode = 'stress';
            this.updateViewModeUI();
            window.app.updateGeometry();
        };

        // Evaluation Mode Toggle
        document.getElementById('mode-fom').onclick = () => {
            window.app.evaluationMode = 'fom';
            this.updateModeUI();
            window.app.updateGeometry();
        };

        document.getElementById('mode-rom').onclick = () => {
            if (window.app.isTrained) {
                window.app.evaluationMode = 'rom';
                this.updateModeUI();
                window.app.updateGeometry();
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
                    console.log("[DOM] Reduction tab active. Chart state:", {
                        appExists: !!window.app,
                        chartsExists: !!(window.app && window.app.charts),
                        energyChartExists: !!(window.app && window.app.charts && window.app.charts.energyChart)
                    });
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

    updateModeUI() {
        const isRom = window.app.evaluationMode === 'rom';
        const fomBtn = document.getElementById('mode-fom');
        const romBtn = document.getElementById('mode-rom');

        if (isRom) {
            romBtn.style.background = 'var(--primary)';
            romBtn.style.color = 'white';
            fomBtn.style.background = 'transparent';
            fomBtn.style.color = '#64748b';
        } else {
            fomBtn.style.background = 'var(--primary)';
            fomBtn.style.color = 'white';
            romBtn.style.background = 'transparent';
            romBtn.style.color = '#64748b';
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
