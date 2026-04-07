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
        const V = [0, 0, 0, 0.5, 0.5, 1, 1, 1]; // Semi-circle (2 arcs)
        
        const n = 9; // Rows
        const m = 5; // Cols
        
        const radius = 15;

        // V-direction Generatrix (Semi-circle in XZ plane)
        const X_gen = [0, radius, radius, radius, 0];
        const Z_gen = [radius, radius, 0, -radius, -radius];
        const V_weights = [1, w, 1, w, 1];

        // U-direction Revolution (Circle in XY plane)
        const u_pts = [
            { x: 1, y: 0 },
            { x: 1, y: 1 },
            { x: 0, y: 1 },
            { x: -1, y: 1 },
            { x: -1, y: 0 },
            { x: -1, y: -1 },
            { x: 0, y: -1 },
            { x: 1, y: -1 },
            { x: 1, y: 0 }
        ];
        const U_weights = [1, w, 1, w, 1, w, 1, w, 1];

        const controlPoints = [];
        const weights = [];
        
        for (let i = 0; i < n; i++) {
            controlPoints[i] = [];
            weights[i] = [];
            
            for (let j = 0; j < m; j++) {
                // Tensor product construction natively produces the correct spherical bounding box
                controlPoints[i][j] = {
                    x: X_gen[j] * u_pts[i].x,
                    y: X_gen[j] * u_pts[i].y,
                    z: Z_gen[j]
                };
                
                weights[i][j] = U_weights[i] * V_weights[j];
            }
        }

        return { p, q, U, V, controlPoints, weights };
    }

    /**
     * Exact NURBS Quadrant for Infinite Plate with Hole Benchmark
     * R: Hole radius, L: Plate dimension
     */
    /**
     * Exact 2-Element NURBS Quadrant (Hughes/Cottrell Benchmark)
     * R: Hole radius, L: Plate dimension
     * This uses a 4x3 control grid with a degenerate corner at (0, L).
     */
    static generatePlateWithHole(R = 1.0, L = 4.0) {
        const p = 2, q = 2;
        const w = Math.cos(Math.PI / 8); // Weight for 22.5-degree quadratic circular sub-arcs
        
        // --- 1. Knot Vectors (2 Elements Angular, 1 Element Radial) ---
        const U = [0, 0, 0, 0.5, 1, 1, 1]; // Angular/Arc
        const V = [0, 0, 0, 1, 1, 1];      // Radial
        
        // --- 2. Control Net (4x3 Grid) ---
        // We define the grid to map a square domain [0, L]x[0, L] to the quarter-plate
        // Inner circle (r=R) at v=0 (j=0)
        // Outer square (x=L or y=L) at v=1 (j=2)
        
        // Point i=0 (0 deg)
        controlPoints[0][0] = { x: R, y: 0, z: 0 };
        controlPoints[0][1] = { x: (R+L)/2, y: 0, z: 0 };
        controlPoints[0][2] = { x: L, y: 0, z: 0 };

        // Point i=1 (22.5 deg)
        controlPoints[1][0] = { x: R, y: R*Math.tan(Math.PI/8), z: 0 };
        controlPoints[1][1] = { x: (R+L)/2, y: (R+L)/2*Math.tan(Math.PI/8), z: 0 };
        controlPoints[1][2] = { x: L, y: L*Math.tan(Math.PI/8), z: 0 };

        // Point i=2 (67.5 deg)
        controlPoints[2][0] = { x: R*Math.tan(Math.PI/8), y: R, z: 0 };
        controlPoints[2][1] = { x: (R+L)/2*Math.tan(Math.PI/8), y: (R+L)/2, z: 0 };
        controlPoints[2][2] = { x: L*Math.tan(Math.PI/8), y: L, z: 0 };

        // Point i=3 (90 deg)
        controlPoints[3][0] = { x: 0, y: R, z: 0 };
        controlPoints[3][1] = { x: 0, y: (R+L)/2, z: 0 };
        controlPoints[3][2] = { x: 0, y: L, z: 0 };

        // --- 3. The "Hughes Sharp Corner" Correction ---
        // To get a sharp corner at (L, L) for a quadratic patch, 
        // the intermediate angular points at the outer radius MUST coalesce.
        controlPoints[1][2] = { x: L, y: L, z: 0 };
        controlPoints[2][2] = { x: L, y: L, z: 0 };

        const weights = [
            [ 1, 1, 1 ],
            [ w, w, w ],
            [ w, w, w ],
            [ 1, 1, 1 ]
        ];

        return { p, q, U, V, controlPoints, weights };
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = NURBSPresets;
} else {
    window.NURBSPresets = NURBSPresets;
}
