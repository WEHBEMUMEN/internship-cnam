/**
 * Phase 5.1a - Geometry Factory
 * Responsible for generating parametric NURBS patches for variable shapes.
 */

class GeometryFactory {
    /**
     * Generates a Plate with a Notch (Quarter Model)
     * @param {number} R  - Notch Radius
     * @param {number} L1 - Plate Width (Horizontal)
     * @param {number} L2 - Plate Height (Vertical)
     */
    static generateNotchPlate(R = 1.0, L1 = 4.0, L2 = 4.0) {
        const p = 2, q = 2;
        const w = Math.cos(Math.PI / 8); // Standard weight for 22.5 deg arcs
        
        // Knot Vectors (2 Elements Angular, 1 Element Radial)
        const U = [0, 0, 0, 0.5, 1, 1, 1]; 
        const V = [0, 0, 0, 1, 1, 1];      

        // 4 Angular Rows (i), 3 Radial Columns (j)
        const controlPoints = [];
        const weights = [];

        // Angles for the quarter notch (0 to 90 degrees)
        const angles = [0, Math.PI/8, 3*Math.PI/8, Math.PI/2];
        const rowWeights = [1, w, w, 1];

        for (let i = 0; i < 4; i++) {
            controlPoints[i] = [];
            weights[i] = [];
            
            const theta = angles[i];
            
            // 1. Inner Point (On the Notch Radius R)
            const P_inner = {
                x: R * Math.cos(theta),
                y: R * Math.sin(theta),
                z: 0
            };

            // 2. Outer Point (On the Plate Boundary L1, L2)
            // We use a mapping that pushes points to the boundary
            let P_outer;
            if (i === 0) { // bottom edge
                P_outer = { x: L1, y: 0, z: 0 };
            } else if (i === 3) { // left edge
                P_outer = { x: 0, y: L2, z: 0 };
            } else if (i === 1) { // 22.5 deg
                P_outer = { x: L1, y: L1 * Math.tan(Math.PI/8), z: 0 };
            } else { // 67.5 deg
                P_outer = { x: L2 * Math.tan(Math.PI/8), y: L2, z: 0 };
            }

            // Adjust for corner degenerate case (The "Hughes" Corner)
            // If it's a square L1=L2, it matches standard benchmarks
            if (i === 1 || i === 2) {
                // For the intermediate angular points, we move towards the "far corner"
                // But we must maintain geometric exactness for the circle.
                // We'll use the user-provided dimensions.
            }

            // Let's refine the outer boundary logic for rectangular plates
            if (i === 1) P_outer = { x: L1, y: L2 * (P_inner.y / R), z: 0 }; // Approx mapping
            if (i === 2) P_outer = { x: L1 * (P_inner.x / R), y: L2, z: 0 };

            // Pure Hughes benchmark logic if L1=L2=L
            if (i === 1) P_outer = { x: L1, y: L1 * Math.tan(Math.PI/8), z: 0 };
            if (i === 2) P_outer = { x: L2 * Math.tan(Math.PI/8), y: L2, z: 0 };

            // 3. Middle Point (Linear blend)
            const P_mid = {
                x: (P_inner.x + P_outer.x) / 2,
                y: (P_inner.y + P_outer.y) / 2,
                z: 0
            };

            controlPoints[i] = [P_inner, P_mid, P_outer];
            weights[i] = [1, rowWeights[i], 1]; // Weights are constant for the radial direction
        }

        // Final Nudge for the Corner (making it a square/rectangle corner)
        // Point [1][2] and [2][2] are the key to the outer square boundary
        // We ensure CP[0][2], [1][2], [2][2], [3][2] form the outer L1 x L2 box
        controlPoints[0][2] = { x: L1, y: 0, z: 0 };
        controlPoints[1][2] = { x: L1, y: L2, z: 0 }; // The corner
        controlPoints[2][2] = { x: L1, y: L2, z: 0 }; // The corner (degenerate)
        controlPoints[3][2] = { x: 0, y: L2, z: 0 };

        // Re-calculate mid-points after corner nudge
        for (let i = 0; i < 4; i++) {
            controlPoints[i][1] = {
                x: (controlPoints[i][0].x + controlPoints[i][2].x) / 2,
                y: (controlPoints[i][0].y + controlPoints[i][2].y) / 2,
                z: 0
            };
        }

        return { 
            p, q, U, V, controlPoints, weights,
            params: { R, L1, L2 }
        };
    }

    /**
     * Applies p and h refinement to a base patch
     */
    static refine(patch, hLevel = 0, pLevel = 0) {
        let refined = JSON.parse(JSON.stringify(patch));
        const engine = new window.NURBS2D();

        // 1. Degree Elevation
        if (pLevel > 0) {
            engine.setDegree(refined, patch.p + pLevel, patch.q + pLevel);
        }

        // 2. Knot Insertion (Global Subdivisions)
        for (let i = 0; i < hLevel; i++) {
            refined = engine.subdivideGlobal(refined);
        }

        return refined;
    }
}

window.GeometryFactory = GeometryFactory;
