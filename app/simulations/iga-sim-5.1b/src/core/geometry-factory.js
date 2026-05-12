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
        const w = Math.cos(Math.PI / 8); 
        
        // Exact logic from Phase 2.2 benchmark NURBSPresets
        const U = [0, 0, 0, 0.5, 1, 1, 1];
        const V = [0, 0, 0, 1, 1, 1];
        
        const controlPoints = [
            [ {x: R, y: 0, z: 0}, {x: (R+L1)/2, y: 0, z: 0}, {x: L1, y: 0, z: 0} ], 
            [ {x: R, y: R*Math.tan(Math.PI/8), z: 0}, {x: L1, y: L1*Math.tan(Math.PI/8), z: 0}, {x: L1, y: L2/2, z: 0} ],
            [ {x: R*Math.tan(Math.PI/8), y: R, z: 0}, {x: L1*Math.tan(Math.PI/8), y: L2, z: 0}, {x: L1/2, y: L2, z: 0} ],
            [ {x: 0, y: R, z: 0}, {x: 0, y: (R+L2)/2, z: 0}, {x: 0, y: L2, z: 0} ]
        ];

        // Degenerate Corner logic from Phase 2.2 benchmark
        controlPoints[1][2] = { x: L1, y: L2, z: 0 };
        controlPoints[2][2] = { x: L1, y: L2, z: 0 };

        // Mid points logic from Phase 2.2 benchmark
        controlPoints[0][1] = { x: (R + L1) / 2, y: 0, z: 0 };
        controlPoints[1][1] = { x: (R + L1) / 2, y: (R * Math.tan(Math.PI/8) + L2) / 2, z: 0 };
        controlPoints[2][1] = { x: (R * Math.tan(Math.PI/8) + L1) / 2, y: (R + L2) / 2, z: 0 };
        controlPoints[3][1] = { x: 0, y: (R + L2) / 2, z: 0 };

        const weights = [
            [ 1, 1, 1 ],
            [ w, w, w ],
            [ w, w, w ],
            [ 1, 1, 1 ]
        ];

        return { 
            p, q, U, V, controlPoints, weights,
            params: { R, L1, L2 },
            elements: GeometryFactory.extractElements({ U, V })
        };
    }

    /**
     * Extracts element definitions (knot spans) for a given patch
     */
    static extractElements(patch) {
        const { U, V } = patch;
        const elements = [];
        for (let i = 0; i < U.length - 1; i++) {
            if (U[i + 1] > U[i]) {
                for (let j = 0; j < V.length - 1; j++) {
                    if (V[j + 1] > V[j]) {
                        elements.push({
                            uRange: [U[i], U[i + 1]],
                            vRange: [V[j], V[j + 1]],
                            spanIdx: [i, j]
                        });
                    }
                }
            }
        }
        return elements;
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

        refined.elements = GeometryFactory.extractElements(refined);
        return refined;
    }
}

window.GeometryFactory = GeometryFactory;
