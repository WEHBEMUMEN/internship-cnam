/**
 * ROM Exporter
 * Phase 2D.1
 */

class ROMExporter {
    static exportSmartBasis(meta, Phi, K1_red, K2_red, Fr) {
        const payload = {
            meta: {
                name: "IGA ROM Basis",
                type: meta.type,
                geometry: meta.geometry,
                analysis: meta.analysis,
                mesh: meta.mesh,
                modesRetained: Phi[0].length
            },
            ROM: {
                Phi: Phi, // Truncated basis matrix
                K1_red: K1_red,
                K2_red: K2_red,
                Fr: Fr // Optional reduced force
            }
        };

        const blob = new Blob([JSON.stringify(payload)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `rom_basis_k${Phi[0].length}.json`;
        a.click();
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { ROMExporter };
} else {
    window.ROMExporter = ROMExporter;
}
