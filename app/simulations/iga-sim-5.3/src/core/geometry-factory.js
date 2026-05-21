/**
 * Phase 5.3 - Geometry Factory
 * Responsible for generating parametric NURBS patches for variable shapes.
 */

class GeometryFactory {
    /**
     * Generates a Cantilever Beam with a Notch at position mu (Phase 5.3)
     * @param {number} L  - Beam Length
     * @param {number} H  - Beam Height
     * @param {number} mu - Notch Center Position
     * @param {number} r  - Notch Depth (from top and bottom)
     */
    static generateNotchedBeam(L = 10, H = 2, mu = 5, r = 0.25) {
        const p = 2, q = 2;
        
        // Use a dense uniform mesh (21 points) for perfect orthogonality
        const nU = 21;
        const nV = 3;
        
        // Uniform open knot vector for U
        const U = [0, 0, 0];
        const numInternalKnots = nU - p - 1; // 21 - 2 - 1 = 18
        for (let i = 1; i <= numInternalKnots; i++) {
            U.push(i / (numInternalKnots + 1));
        }
        U.push(1, 1, 1);
        
        // 3 points in V (Height)
        const V = [0, 0, 0, 1, 1, 1];
        
        const controlPoints = [];
        const weights = [];
        
        for (let i = 0; i < nU; i++) {
            controlPoints[i] = [];
            weights[i] = [];
            
            // Uniform X spacing
            let x = (i / (nU - 1)) * L;
            
            for (let j = 0; j < nV; j++) {
                // Centered at Y=0
                let y = (j === 0) ? -H/2 : (j === 1 ? 0 : H/2);
                
                // Defect shape: Smooth Gaussian notch centered at mu
                let dist = Math.abs(x - mu);
                // The notch width is roughly proportional to its depth r
                let sigma = Math.max(r, 0.4); 
                
                if (dist < sigma * 3) {
                    // Calculate displacement
                    let disp = r * Math.exp(-0.5 * Math.pow(dist / sigma, 2));
                    
                    if (j === 0) y += disp; // Bottom notch goes up
                    if (j === 2) y -= disp; // Top notch goes down
                }
                
                controlPoints[i][j] = { x, y, z: 0 };
                weights[i][j] = 1;
            }
        }
        
        return {
            p, q, U, V, controlPoints, weights,
            params: { L, H, mu, r },
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
