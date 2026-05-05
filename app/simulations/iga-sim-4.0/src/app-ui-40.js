/**
 * Phase 4.0 UI Logic
 */

document.addEventListener('DOMContentLoaded', () => {
    const app = window.app;

    // Range Inputs
    const inputs = [
        { id: 'input-time', valId: 'time-val', suffix: 's' },
        { id: 'input-dt', valId: 'dt-val', suffix: 's' },
        { id: 'input-alpha', valId: 'alpha-val', suffix: '' },
        { id: 'input-beta', valId: 'beta-val', suffix: '' },
        { id: 'input-f0', valId: 'f0-val', suffix: 'N' }
    ];

    inputs.forEach(cfg => {
        const el = document.getElementById(cfg.id);
        const valEl = document.getElementById(cfg.valId);
        el.addEventListener('input', () => {
            valEl.innerText = el.value + cfg.suffix;
        });
    });

    // Run Button
    const btnRun = document.getElementById('btn-run');
    btnRun.addEventListener('click', async () => {
        if (app.isRunning) {
            app.isRunning = false;
            btnRun.innerHTML = '<i class="fa-solid fa-play"></i> Run Transient FOM';
            document.getElementById('run-label').innerText = 'Stopped';
            document.getElementById('run-dot').style.background = '#f43f5e';
        } else {
            btnRun.innerHTML = '<i class="fa-solid fa-stop"></i> Stop Simulation';
            document.getElementById('run-label').innerText = 'Running';
            document.getElementById('run-dot').style.background = '#10b981';
            await app.runSimulation();
            btnRun.innerHTML = '<i class="fa-solid fa-play"></i> Run Transient FOM';
        }
    });
});
