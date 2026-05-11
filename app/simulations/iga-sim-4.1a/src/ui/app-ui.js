/**
 * Phase 4.1 UI Logic
 */

document.addEventListener('DOMContentLoaded', () => {
    const app = window.app;

    // Range Inputs
    const inputs = [
        { id: 'input-time', valId: 'time-val', suffix: 's' },
        { id: 'input-steps', valId: 'steps-val', suffix: '' },
        { id: 'input-defscale', valId: 'defscale-val', suffix: 'x' },
        { id: 'input-k', valId: 'k-val', suffix: '' },
        { id: 'input-m', valId: 'm-val', suffix: '' },
        { id: 'input-h', valId: 'h-val', suffix: '' },
        { id: 'input-p', valId: 'p-val', suffix: '' },
        { id: 'input-fmin', valId: 'fmin-val', suffix: 'N' },
        { id: 'input-fmax', valId: 'fmax-val', suffix: 'N' },
        { id: 'input-samples', valId: 'samples-val', suffix: '' }
    ];

    inputs.forEach(cfg => {
        const el = document.getElementById(cfg.id);
        const valEl = document.getElementById(cfg.valId);
        if (!el) return;
        el.addEventListener('input', () => {
            let displayVal = el.value;
            if (cfg.id === 'input-m') {
                displayVal = "1e-" + el.value;
                if (!app.isTraining && app.snapshots.length > 0) {
                    app.recomputeECSW(parseInt(el.value));
                }
            }
            valEl.innerText = displayVal + cfg.suffix;

            if (cfg.id === 'input-defscale') {
                app.viz.defScale = parseFloat(el.value);
                if (app.dyn && app.dyn.u) app.viz.updateMesh(app.dyn.u);
            } else if (cfg.id === 'input-h' || cfg.id === 'input-p') {
                const hVal = parseInt(document.getElementById('input-h').value);
                const pVal = parseInt(document.getElementById('input-p').value);
                
                // Show warning if high intensity
                const warning = document.getElementById('mesh-warning');
                if (warning) warning.style.display = hVal >= 4 ? 'block' : 'none';
                
                app.updateRefinement(hVal, pVal);
            }
        });
    });

    // Control Points Toggle
    const checkCP = document.getElementById('check-cp');
    const toggleTrigger = document.getElementById('toggle-cp-trigger');
    const toggleDot = document.getElementById('toggle-dot');
    const toggleBg = document.getElementById('toggle-bg');
    
    if (toggleTrigger) {
        toggleTrigger.addEventListener('click', () => {
            checkCP.checked = !checkCP.checked;
            app.viz.showCP = checkCP.checked;
            toggleDot.style.transform = checkCP.checked ? 'translateX(20px)' : 'translateX(0)';
            toggleBg.style.background = checkCP.checked ? 'var(--primary)' : '#334155';
            if (app.dyn && app.dyn.u) app.viz.updateMesh(app.dyn.u);
        });
    }

    // Toggle Active Elements
    const toggleActive = document.getElementById('toggle-active');
    if (toggleActive) {
        toggleActive.addEventListener('change', (e) => {
            app.viz.showActiveElements = e.target.checked;
            if (app.dyn && app.dyn.u) app.viz.updateMesh(app.dyn.u);
        });
    }

    // Draggable Logic (from 4.0)
    let activePanel = null;
    let offset = [0, 0];
    let isDown = false;

    document.querySelectorAll('.glass-panel').forEach(panel => {
        panel.addEventListener('mousedown', (e) => {
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'BUTTON') return;
            isDown = true;
            activePanel = panel;
            offset = [panel.offsetLeft - e.clientX, panel.offsetTop - e.clientY];
            document.querySelectorAll('.glass-panel').forEach(p => p.style.zIndex = "10");
            panel.style.zIndex = "100";
        });
    });

    document.addEventListener('mouseup', () => { isDown = false; });
    document.addEventListener('mousemove', (e) => {
        if (isDown && activePanel) {
            activePanel.style.left = (e.clientX + offset[0]) + 'px';
            activePanel.style.top = (e.clientY + offset[1]) + 'px';
            activePanel.style.right = 'auto';
            activePanel.style.bottom = 'auto';
        }
    });

    // Run Training
    document.getElementById('btn-train').addEventListener('click', async () => {
        await app.runTraining();
        document.getElementById('input-k').disabled = false;
        document.getElementById('input-m').disabled = false;
    });

    const btnAudit = document.getElementById('btn-audit');
    if (btnAudit) {
        btnAudit.addEventListener('click', () => {
            window.Audit41a.run(window.app);
        });
    }

    // Export Package
    const btnExport = document.getElementById('btn-export');
    if (btnExport) {
        btnExport.addEventListener('click', () => {
            app.trainer.exportPackage();
        });
    }

    // Verify Package
    const btnVerify = document.getElementById('btn-verify');
    if (btnVerify) {
        btnVerify.addEventListener('click', () => {
            const report = app.trainer.verifyPackage();
            showVerificationModal(report);
        });
    }

    function showVerificationModal(report) {
        const overlay = document.createElement('div');
        overlay.style = "position:fixed; inset:0; background:rgba(15,23,42,0.8); backdrop-filter:blur(8px); display:flex; align-items:center; justify-content:center; z-index:2000;";
        
        const modal = document.createElement('div');
        modal.className = "glass-panel";
        modal.style = "width:420px; padding:30px; cursor:default; position:relative; border-color:var(--primary);";
        
        let html = `
            <div style="font-size:0.7rem; font-weight:800; text-transform:uppercase; letter-spacing:2px; color:var(--primary); margin-bottom:20px; text-align:center;">ROM Verification Audit</div>
            <div style="display:flex; flex-direction:column; gap:12px; margin-bottom:25px;">
        `;

        report.checks.forEach(c => {
            const color = c.status === "OK" ? "#10b981" : (c.status === "FAIL" ? "#f43f5e" : "#f59e0b");
            const icon = c.status === "OK" ? "fa-circle-check" : (c.status === "FAIL" ? "fa-circle-xmark" : "fa-triangle-exclamation");
            
            html += `
                <div style="background:rgba(30,41,59,0.5); border-radius:12px; padding:12px; border:1px solid rgba(255,255,255,0.05); display:flex; justify-content:space-between; align-items:center;">
                    <div>
                        <div style="font-size:0.65rem; color:#94a3b8; font-weight:700; text-transform:uppercase; letter-spacing:1px;">${c.name}</div>
                        <div style="font-size:0.8rem; font-weight:600; color:#f1f5f9; margin-top:2px;">${c.detail}</div>
                    </div>
                    <div style="color:${color}; font-size:1.1rem;"><i class="fa-solid ${icon}"></i></div>
                </div>
            `;
        });

        html += `
            </div>
            <button id="close-modal" class="btn-action" style="background:#1e293b; border:1px solid var(--border);">Close Audit Report</button>
        `;

        modal.innerHTML = html;
        overlay.appendChild(modal);
        document.body.appendChild(overlay);

        modal.querySelector('#close-modal').onclick = () => overlay.remove();
        overlay.onclick = (e) => { if(e.target === overlay) overlay.remove(); };
    }
});
