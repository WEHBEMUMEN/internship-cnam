/**
 * NURBS Presets for Phase 2.2
 * Contains specialized geometries like exact Spheres and Toruses.
 */

class NURBSPresets {
    static generateSheet() {
        const patch = {
            p: 2, q: 2,
            U: [0, 0, 0, 0.5, 1, 1, 1],
            V: [0, 0, 0, 0.5, 1, 1, 1],
            weights: Array(4).fill().map(() => Array(4).fill(1)),
            controlPoints: []
        };

        for (let i = 0; i < 4; i++) {
            patch.controlPoints[i] = [];
            for (let j = 0; j < 4; j++) {
                patch.controlPoints[i][j] = {
                    x: (i - 1.5) * 10,
                    y: (j - 1.5) * 10,
                    z: 0
                };
            }
        }
        
        // Default deformation
        patch.controlPoints[0][0].z = 5; patch.controlPoints[0][3].z = -5;
        patch.controlPoints[3][0].z = -5; patch.controlPoints[3][3].z = 5;
        patch.controlPoints[1][1].z = 15; patch.controlPoints[2][2].z = 15;
        
        return patch;
    }

    /**
     * Exact NURBS Sphere - Single-Patch Surface of Revolution
     * Constructed by revolving a semi-circular NURBS curve (degree p=2) around a central axis.
     * Note: This method generates "poles" (singularities) at the top and bottom.
     */
    static generateSphere() {
        const p = 2; // Quadratic, minimum required for circular exactness
        const q = 2;
        
        // Exact Circle/Revolve logic requires rational weights.
        // For a 90-degree circular arc (degree 2), weights must be: { 1, sqrt(2)/2, 1 }
        const w = Math.sqrt(2) / 2;
        
        // Open (clamped) knot vectors: First and last knots repeated p+1 times
        const U = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]; // Full revolution (4 arcs)
        const V = [0, 0, 0, 0.5, 1, 1, 1]; // Semi-circle (2 arcs)
        
        const n = 9; // Rows
        const m = 3; // Cols
        
        const controlPoints = [];
        const weights = [];
        
        const radius = 15;
        
        for (let i = 0; i < n; i++) {
            controlPoints[i] = [];
            weights[i] = [];
            
            // Angle in U (Full circle, 360 deg)
            const phi = (i * Math.PI) / 4; 
            const weightU = (i % 2 === 1) ? w : 1;
            
            for (let j = 0; j < m; j++) {
                // Angle in V (Semi circle, 180 deg)
                const theta = (j * Math.PI) / 2;
                const weightV = (j % 2 === 1) ? w : 1;
                
                // Sphere parametric eq:
                // x = R * cos(phi) * sin(theta)
                // y = R * sin(phi) * sin(theta)
                // z = R * cos(theta)
                
                // Adjusting for the NURBS "Pull" of the middle points
                let R_eff = radius;
                if (i % 2 === 1) R_eff /= weightU;
                if (j === 1) R_eff /= weightV;

                controlPoints[i][j] = {
                    x: R_eff * Math.cos(phi) * Math.sin(theta),
                    y: R_eff * Math.sin(phi) * Math.sin(theta),
                    z: R_eff * Math.cos(theta)
                };
                
                weights[i][j] = weightU * weightV;
            }
        }

        return { p, q, U, V, controlPoints, weights };
    }
}

window.NURBSPresets = NURBSPresets;
