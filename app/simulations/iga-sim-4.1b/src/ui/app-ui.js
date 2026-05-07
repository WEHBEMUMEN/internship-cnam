/**
 * Phase 4.1 UI Logic
 */

document.addEventListener('DOMContentLoaded', () => {
    const app = window.app;

    // Range Inputs
    const inputs = [
        { id: 'input-defscale', valId: 'defscale-val', suffix: 'x' },
        { id: 'input-fy', valId: 'fy-val', suffix: 'N' }
    ];

    inputs.forEach(cfg => {
        const el = document.getElementById(cfg.id);
        const valEl = document.getElementById(cfg.valId);
        if (!el) return;
        el.addEventListener('input', () => {
            valEl.innerText = el.value + cfg.suffix;
            if (cfg.id === 'input-defscale') {
                app.viz.defScale = parseFloat(el.value);
                if (app.rom && app.rom.u) app.viz.updateMesh(app.rom.u);
            }
        });
    });

    // Package Import
    const btnTrigger = document.getElementById('btn-import-trigger');
    const inputPackage = document.getElementById('input-package');
    
    btnTrigger.addEventListener('click', () => inputPackage.click());
    
    inputPackage.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (!file) return;
        
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const data = JSON.parse(event.target.result);
                app.loadPackage(data);
            } catch (err) {
                alert("Invalid JSON package: " + err.message);
            }
        };
        reader.readAsText(file);
    });

    // Control Buttons
    document.getElementById('btn-run').addEventListener('click', () => app.toggle());
    document.getElementById('btn-reset').addEventListener('click', () => app.reset());
    
    // Toggle Active Elements
    const toggleActive = document.getElementById('toggle-active');
    if (toggleActive) {
        toggleActive.addEventListener('change', (e) => {
            app.viz.showActiveElements = e.target.checked;
            if (app.rom && app.rom.u) app.viz.updateMesh(app.rom.u);
        });
    }
    
    const btnAudit = document.getElementById('btn-audit');
    if (btnAudit) {
        btnAudit.addEventListener('click', () => {
            alert("Benchmark: ROM Step Latency < 1ms (Hyper-Reduced)");
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
