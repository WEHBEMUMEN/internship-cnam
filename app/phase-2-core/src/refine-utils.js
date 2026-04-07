/**
 * Refinement Utilities for IGA
 * Provides a standardized way to apply h, p, and k-refinement stacks.
 */
class RefineUtils {
    /**
     * Apply a refinement stack (degree followed by subdivision)
     * @param {NURBS2D} engine 
     * @param {Object} patch 
     * @param {Object} options { h, p } 
     */
    static apply(engine, patch, options = { h: 0, p: 2 }) {
        const { h, p } = options;
        
        // 1. Polynomial Elevation (p-refinement)
        // Note: engine.elevateDirection mutates in-place
        if (p > patch.p) {
            const deltaP = p - patch.p;
            for (let i = 0; i < deltaP; i++) {
                engine.elevateDirection(patch, 'U');
                engine.elevateDirection(patch, 'V');
            }
        }
        
        // 2. Global Subdivision (h-refinement)
        // Note: engine.subdivideGlobal mutates in-place (and returns patch)
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

if (typeof module !== 'undefined' && module.exports) {
    module.exports = RefineUtils;
} else {
    window.RefineUtils = RefineUtils;
}
