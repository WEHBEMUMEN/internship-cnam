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
    static generatePlateWithHole(R = 1.0, L = 4.0) {
        const p = 2, q = 2;
        const w = 1 / Math.sqrt(2); 
        
        // --- 1. Geometry Setup (Top-Right Quadrant) ---
        // U direction: Radial (Hole to Outer Edge)
        // V direction: Angular (0 to 90 degrees)
        const U = [0, 0, 0, 1, 1, 1];
        const V = [0, 0, 0, 1, 1, 1];
        
        const controlPoints = [];
        const weights = [];

        // We construct 3 layers of 3 points (Quadratic 3x3 Grid)
        // Layer 0: Hole radius R
        // Layer 1: Intermediate
        // Layer 2: Outer edge L
        
        // Inner Circular Boundary (u=0)
        const inner = [
            { x: R, y: 0, z: 0 },
            { x: R, y: R, z: 0 }, // Weight w
            { x: 0, y: R, z: 0 }
        ];

        // Outer Square Boundary (u=1)
        const outer = [
            { x: L, y: 0, z: 0 },
            { x: L, y: L, z: 0 },
            { x: 0, y: L, z: 0 }
        ];

        // Mid Layer (Linear interpolate R and L)
        const mid = inner.map((p, i) => ({
            x: p.x + (outer[i].x - p.x) * 0.5,
            y: p.y + (outer[i].y - p.y) * 0.5,
            z: 0
        }));

        const grid = [inner, mid, outer];
        for(let i=0; i<3; i++) {
            controlPoints[i] = grid[i];
            weights[i] = [1, w, 1];
        }

        const patch = { p, q, U, V, controlPoints, weights };

        // --- 2. Initial p-Refinement (Ensure sufficient nodes for the plot) ---
        // Increase resolution to p=3, q=3 for smoother plotting and assembly
        const engine = new NURBS2D();
        engine.elevateDirection(patch, 'U');
        engine.elevateDirection(patch, 'V');

        return patch;
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = NURBSPresets;
} else {
    window.NURBSPresets = NURBSPresets;
}
