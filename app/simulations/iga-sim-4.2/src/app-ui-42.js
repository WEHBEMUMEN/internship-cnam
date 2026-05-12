/**
 * Phase 4.0 UI Logic
 */

document.addEventListener('DOMContentLoaded', () => {
    const app = window.app;

    // Range Inputs
    const inputs = [
        { id: 'input-dt', valId: 'dt-val', suffix: 's' },
        { id: 'input-alpha', valId: 'alpha-val', suffix: '' },
        { id: 'input-beta', valId: 'beta-val', suffix: '' },
        { id: 'input-f0', valId: 'f0-val', suffix: 'N' },
        { id: 'input-res', valId: 'res-val', suffix: 'x' },
        { id: 'input-h', valId: 'h-val', suffix: '' },
        { id: 'input-p', valId: 'p-val', suffix: '' },
        { id: 'input-timescale', valId: 'timescale-val', suffix: 'x' },
        { id: 'input-defscale', valId: 'defscale-val', suffix: 'x' },
        { id: 'input-young', valId: 'young-val', suffix: '' },
        { id: 'input-poisson', valId: 'poisson-val', suffix: '' },
        { id: 'input-rho', valId: 'rho-val', suffix: 'kg/m³' },
        { id: 'input-gravity', valId: 'gravity-val', suffix: ' m/s²' }
    ];

    inputs.forEach(cfg => {
        const el = document.getElementById(cfg.id);
        const valEl = document.getElementById(cfg.valId);
        if (!el) return;
        el.addEventListener('input', () => {
            if (cfg.id === 'input-res') {
                valEl.innerText = `${el.value}x${el.value}`;
                app.viz.meshRes = parseInt(el.value);
                if (!app.isRunning) app.viz.updateMesh(app.dyn.u);
            } else if (cfg.id === 'input-h' || cfg.id === 'input-p') {
                valEl.innerText = el.value;
                const hVal = parseInt(document.getElementById('input-h').value);
                const pVal = parseInt(document.getElementById('input-p').value);
                app.updateRefinement(hVal, pVal);
                // Update Run button state
                const btnRun = document.getElementById('btn-run');
                if (btnRun) btnRun.innerHTML = '<i class="fa-solid fa-play"></i> Play';
                const runLabel = document.getElementById('run-label');
                if (runLabel) runLabel.innerText = 'Ready';
                const runDot = document.getElementById('run-dot');
                if (runDot) runDot.style.background = '#10b981';
            } else if (cfg.id === 'input-alpha') {
                valEl.innerText = el.value;
                app.dyn.alpha = parseFloat(el.value);
            } else if (cfg.id === 'input-beta') {
                valEl.innerText = el.value;
                app.dyn.betaR = parseFloat(el.value);
            } else if (cfg.id === 'input-defscale') {
                valEl.innerText = el.value + 'x';
                app.viz.defScale = parseFloat(el.value);
                if (!app.isRunning) app.viz.updateMesh(app.dyn.u);
            } else if (cfg.id === 'input-young') {
                const val = parseFloat(el.value);
                valEl.innerText = val.toExponential(2);
                app.E = val;
                app.fom.E = val;
            } else if (cfg.id === 'input-poisson') {
                const val = parseFloat(el.value);
                valEl.innerText = val.toFixed(2);
                app.nu = val;
                app.fom.nu = val;
            } else if (cfg.id === 'input-rho') {
                const val = parseFloat(el.value);
                valEl.innerText = el.value + 'kg/m³';
                app.rho = val;
                app.dyn.rho = val;
                app.dyn.assembleMass();
            } else {
                valEl.innerText = el.value + cfg.suffix;
            }
        });
    });

    // Control Points Toggle
    const checkCP = document.getElementById('check-cp');
    const toggleTrigger = document.getElementById('toggle-cp-trigger');
    const toggleDot = document.getElementById('toggle-dot');
    const toggleBg = document.getElementById('toggle-bg');
    
    toggleTrigger.addEventListener('click', () => {
        checkCP.checked = !checkCP.checked;
        app.viz.showCP = checkCP.checked;
        toggleDot.style.transform = checkCP.checked ? 'translateX(20px)' : 'translateX(0)';
        toggleBg.style.background = checkCP.checked ? 'var(--primary)' : '#334155';
        if (!app.isRunning) app.viz.updateMesh(app.dyn.u);
    });

    // Draggable Panels
    const panels = document.querySelectorAll('.glass-panel');
    panels.forEach(panel => {
        let isDragging = false;
        let startX, startY, initialX, initialY;

        panel.addEventListener('mousedown', (e) => {
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.tagName === 'BUTTON' || e.target.closest('.toggle-container')) return;
            isDragging = true;
            startX = e.clientX;
            startY = e.clientY;
            initialX = panel.offsetLeft;
            initialY = panel.offsetTop;
            panel.style.zIndex = 100;
            // Bring to front
            panels.forEach(p => { if (p !== panel) p.style.zIndex = 10; });
        });

        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const dx = e.clientX - startX;
            const dy = e.clientY - startY;
            panel.style.left = `${initialX + dx}px`;
            panel.style.top = `${initialY + dy}px`;
            panel.style.right = 'auto'; // Break alignment
            panel.style.bottom = 'auto'; 
        });

        document.addEventListener('mouseup', () => {
            isDragging = false;
        });
    });

    // Pause/Play Toggle Button
    const btnRun = document.getElementById('btn-run');
    btnRun.addEventListener('click', async () => {
        if (app.isRunning) {
            app.isRunning = false;
            btnRun.innerHTML = '<i class="fa-solid fa-play"></i> Play';
            document.getElementById('run-label').innerText = 'Paused';
            document.getElementById('run-dot').style.background = '#f59e0b';
        } else {
            btnRun.innerHTML = '<i class="fa-solid fa-pause"></i> Pause';
            document.getElementById('run-label').innerText = 'Running';
            document.getElementById('run-dot').style.background = '#10b981';
            await app.runSimulation();
        }
    });

    // Reset Button
    const btnReset = document.getElementById('btn-reset');
    btnReset.addEventListener('click', () => {
        app.isRunning = false;
        app.currentTime = 0;
        app.dyn.u.fill(0);
        app.dyn.v.fill(0);
        app.dyn.a.fill(0);
        app.trace = [];
        app.mainChart.data.labels = [];
        app.mainChart.data.datasets[0].data = [];
        app.mainChart.update();
        app.responseChart.data.labels = [];
        app.responseChart.data.datasets[0].data = [];
        app.responseChart.update();
        app.viz.updateMesh(app.dyn.u);
        
        btnRun.innerHTML = '<i class="fa-solid fa-play"></i> Play';
        document.getElementById('run-label').innerText = 'Ready';
        document.getElementById('run-dot').style.background = '#10b981';
    });

    // CSM Case Selector
    const caseSelect = document.getElementById('input-csm-case');
    caseSelect.addEventListener('change', () => {
        const csmCase = caseSelect.value;
        const youngInput = document.getElementById('input-young');
        const gravityInput = document.getElementById('input-gravity');
        const targetLabel = document.getElementById('target-y');
        const targetFreqLabel = document.getElementById('target-f');
        
        if (csmCase === 'csm1') {
            app.targetY = -7.187;
            youngInput.value = 1400000;
            gravityInput.value = 2.0;
            targetLabel.innerText = '-7.187 mm';
            targetFreqLabel.innerText = '0 Hz (Static)';
        } else if (csmCase === 'csm2') {
            app.targetY = -1.897;
            youngInput.value = 5600000;
            gravityInput.value = 2.0;
            targetLabel.innerText = '-1.897 mm';
            targetFreqLabel.innerText = '0 Hz (Static)';
        } else if (csmCase === 'csm3') {
            app.targetY = -3.421; // Transient mean
            youngInput.value = 1400000;
            gravityInput.value = 2.0;
            targetLabel.innerText = 'Mean: -3.42 mm';
            targetFreqLabel.innerText = '1.0995 Hz';
        }
        
        // Trigger inputs
        youngInput.dispatchEvent(new Event('input'));
        gravityInput.dispatchEvent(new Event('input'));
        
        if (!app.isRunning) app.viz.updateMesh(app.dyn.u);
    });

    // Audit Button
    const btnAudit = document.getElementById('btn-audit');
    if (btnAudit) {
        btnAudit.addEventListener('click', () => {
            window.Audit42.run(window.app);
        });
    }
});
