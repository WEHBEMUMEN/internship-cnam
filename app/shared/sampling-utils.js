/**
 * Sampling Utilities for IGA-ROM
 * Phase 2C | Offline State Generation
 */

class SamplingUtils {
    /**
     * Latin Hypercube Sampling (LHS)
     * @param {Array} dimensions Array of { name, min, max }
     * @param {number} nSamples Number of samples to generate
     * @returns {Array} Array of sample objects
     */
    static generateLHS(dimensions, nSamples) {
        const nDim = dimensions.length;
        const result = [];

        // 1. Create bins for each dimension
        const bins = dimensions.map(dim => {
            const stride = (dim.max - dim.min) / nSamples;
            const b = [];
            for (let i = 0; i < nSamples; i++) {
                // Random value within the bin
                const val = dim.min + (i + Math.random()) * stride;
                b.push(val);
            }
            // 2. Shuffle each dimension separately
            return this.shuffle(b);
        });

        // 3. Combine bins into samples
        for (let i = 0; i < nSamples; i++) {
            const sample = {};
            dimensions.forEach((dim, dIdx) => {
                sample[dim.name] = bins[dIdx][i];
            });
            result.push(sample);
        }

        return result;
    }

    /**
     * Fisher-Yates Shuffle
     */
    static shuffle(array) {
        for (let i = array.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [array[i], array[array[j] ? j : i]] = [array[j], array[i]];
        }
        return array;
    }
}

// Export for browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { SamplingUtils };
} else {
    window.SamplingUtils = SamplingUtils;
}
