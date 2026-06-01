/**
 * Refinement Utilities for IGA shared library.
 * Centralizes degree elevations (p-refinement), global subdivisions (h-refinement),
 * and k-refinement stacks.
 */

class RefineUtils {
    /**
     * Apply a refinement stack (degree followed by subdivision)
     * @param {NURBS2D|NURBSCore} engine The NURBS engine containing subdivision/elevation code
     * @param {Object} patch NURBS patch object
     * @param {Object} options { h, p } parameters
     * @returns {Object} Refined patch
     */
    static apply(engine, patch, options = { h: 0, p: 2 }) {
        const { h, p } = options;
        
        // 1. Polynomial Elevation (p-refinement)
        if (p > patch.p) {
            const deltaP = p - patch.p;
            for (let i = 0; i < deltaP; i++) {
                engine.elevateDirection(patch, 'U');
                engine.elevateDirection(patch, 'V');
            }
        }
        
        // 2. Global Subdivision (h-refinement)
        if (h > 0) {
            for (let i = 0; i < h; i++) {
                engine.subdivideGlobal(patch);
            }
        }
        
        return patch;
    }

    /**
     * K-Refinement: Standard Elevate-then-Insert Strategy
     * Ensures C^{p-1} continuity throughout.
     */
    static applyKRef(engine, patch, { hDelta, pDelta }) {
        // First elevate degree (p-ref)
        for (let i = 0; i < pDelta; i++) {
            engine.elevateDirection(patch, 'U');
            engine.elevateDirection(patch, 'V');
        }
        // Then subdivide (h-ref)
        for (let i = 0; i < hDelta; i++) {
            engine.subdivideGlobal(patch);
        }
        return patch;
    }
}

// Export for environment compatibility
if (typeof module !== 'undefined' && module.exports) {
    module.exports = RefineUtils;
} else {
    window.RefineUtils = RefineUtils;
}
