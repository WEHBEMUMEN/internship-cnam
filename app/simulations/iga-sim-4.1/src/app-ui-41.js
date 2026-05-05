/**
 * Phase 4.1 UI Logic
 */

document.addEventListener('DOMContentLoaded', () => {
    const app = window.app;

    // Range Inputs
    const inputs = [
        { id: 'input-time', valId: 'time-val', suffix: 's' },
        { id: 'input-steps', valId: 'steps-val', suffix: '' },
        { id: 'input-k', valId: 'k-val', suffix: '' },
        { id: 'input-m', valId: 'm-val', suffix: '' }
    ];

    inputs.forEach(cfg => {
        const el = document.getElementById(cfg.id);
        const valEl = document.getElementById(cfg.valId);
        if (!el) return;
        el.addEventListener('input', () => {
            let displayVal = el.value;
            if (cfg.id === 'input-m') displayVal = "1e-" + el.value;
            valEl.innerText = displayVal + cfg.suffix;
        });
    });

    // Run Training
    document.getElementById('btn-train').addEventListener('click', async () => {
        await app.runTraining();
        document.getElementById('input-k').disabled = false;
        document.getElementById('input-m').disabled = false;
    });

    // Export Package
    document.getElementById('btn-export').addEventListener('click', () => {
        app.trainer.exportPackage();
    });

    // Verify Package
    document.getElementById('btn-verify').addEventListener('click', () => {
        const isValid = app.trainer.verifyPackage();
        if (isValid) {
            alert("ROM Package Verification Successful!\nContents are consistent and ready for online use.");
        } else {
            alert("ROM Package Verification Failed. Check console for details.");
        }
    });
});
