/**
 * Phase 5.1a - Visuals Engine (Three.js)
 */

class VisualsEngine {
    constructor(containerId, nurbs, solver) {
        this.container = document.getElementById(containerId);
        this.nurbs = nurbs;
        this.solver = solver;
        this.defScale = 50.0;
        this.showControlNet = false;
        this.showBCs = true; 
        this.viewMode = 'displacement'; 
        
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x020617); // Exact 2B.2 dark
        
        this.camera = new THREE.PerspectiveCamera(45, this.container.clientWidth / this.container.clientHeight, 0.1, 1000);
        this.camera.position.set(2.5, 2.5, 12); // Center of 4x4 plate
        this.camera.lookAt(2.5, 2.5, 0);

        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.screenSpacePanning = true;
        this.controls.mouseButtons = {
            LEFT: THREE.MOUSE.PAN,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.PAN
        };
        this.controls.target.set(2.5, 2.5, 0);
        this.controls.update();

        // Lighting
        const ambient = new THREE.AmbientLight(0xffffff, 0.8);
        this.scene.add(ambient);
        const directional = new THREE.DirectionalLight(0xffffff, 0.8);
        directional.position.set(10, 10, 20);
        this.scene.add(directional);

        // Objects
        this.mesh = null;
        this.wireframe = null;
        this.cpNet = null;

        window.addEventListener('resize', () => this.onResize());
        this.animate();
    }

    onResize() {
        this.camera.aspect = this.container.clientWidth / this.container.clientHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
    }

    /**
     * Updates the 3D representation of the NURBS patch
     */
    update(patch, displacement = null) {
        if (this.mesh) this.scene.remove(this.mesh);
        if (this.wireframe) this.scene.remove(this.wireframe);
        if (this.cpNet) this.scene.remove(this.cpNet);

        const segments = 48; // Higher resolution for stress
        const geometry = new THREE.PlaneGeometry(1, 1, segments, segments);
        const positions = geometry.attributes.position;
        const colors = new Float32Array(positions.count * 3);
        
        let maxVal = 0;
        const vals = new Float32Array(positions.count);

        // 1. Evaluate Surface and Scalar Field
        console.log("[Visuals] Starting update. Mode:", this.viewMode, "Solver:", !!this.solver);
        
        // Pre-calculate max for normalization (matching 2B.2 logic)
        if (displacement) {
            if (this.viewMode === 'displacement') {
                for (let i = 0; i < displacement.length; i += 2) {
                    // Match 2B.2: maxVal for the color bar is the Resultant Magnitude
                    const mag = Math.sqrt(displacement[i]**2 + displacement[i+1]**2);
                    if (mag > maxVal) maxVal = mag;
                }
            } else if (this.viewMode === 'stress' && this.solver) {
                // Sample stress capped at v=0.95 to avoid degenerate corner artifacts (2B.2 logic)
                const res = 15;
                for (let i = 0; i <= res; i++) {
                    for (let j = 0; j <= res; j++) {
                        const u = i / res;
                        const v = (j / res) * 0.95;
                        const s = this.solver.getNumericalStress(patch, displacement, u, v, this.solver.E, this.solver.nu);
                        if (s.vonMises > maxVal) maxVal = s.vonMises;
                    }
                }
            }
        }
        if (maxVal === 0) maxVal = 1.0;

        for (let i = 0; i < positions.count; i++) {
            const u = geometry.attributes.uv.getX(i);
            const v = geometry.attributes.uv.getY(i);
            
            const pos = this.nurbs.evaluateSurface(patch, u, v);
            
            if (displacement) {
                const u_disp = this.evaluateDisplacement(patch, u, v, displacement);
                pos.x += u_disp.x * this.defScale;
                pos.y += u_disp.y * this.defScale;

                if (this.viewMode === 'stress' && this.solver) {
                    try {
                        const s = this.solver.getNumericalStress(patch, displacement, u, v, this.solver.E, this.solver.nu);
                        vals[i] = s.vonMises;
                    } catch (err) {
                        if (i === 0) console.error("[Visuals] Stress Eval Error:", err);
                    }
                } else if (this.viewMode === 'displacement') {
                    // Match 2B.2: Use absolute X-displacement for coloring
                    vals[i] = Math.abs(u_disp.x);
                }
            }
            positions.setXYZ(i, pos.x, pos.y, pos.z);
        }

        // 2. Map Colors (Exact 2B.2 Heatmap Stops)
        for (let i = 0; i < positions.count; i++) {
            const color = this.getJetColor(vals[i], 0, maxVal);
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
        }
        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color: 0xffffff,
            vertexColors: true,
            side: THREE.DoubleSide,
            flatShading: false,
            transparent: true,
            opacity: 0.9,
            shininess: 100
        });

        this.mesh = new THREE.Mesh(geometry, material);
        this.scene.add(this.mesh);

        const wireMaterial = new THREE.MeshBasicMaterial({ color: 0x475569, wireframe: true, transparent: true, opacity: 0.3 });
        this.wireframe = new THREE.Mesh(geometry, wireMaterial);
        this.scene.add(this.wireframe);

        // 3. ECSW Active Mesh
        if (window.app && window.app.ecswData) {
            this.renderECSWActiveMesh(patch, window.app.ecswData);
        }

        this.renderControlNet(patch);
        
        // Update Color Bar
        const title = this.viewMode === 'stress' ? 'Von Mises (MPa)' : 'Disp-X (mm)';
        this.updateColorBar(0, maxVal, title);
    }

    updateColorBar(min, max, title) {
        const container = document.getElementById('colorbar-container');
        if (!container) return;
        
        if (max === 0) {
            container.style.display = 'none';
            return;
        }
        
        container.style.display = 'flex';
        document.getElementById('cb-max').textContent = max.toFixed(3);
        document.getElementById('cb-mid').textContent = ((max + min) / 2).toFixed(3);
        document.getElementById('cb-min').textContent = min.toFixed(3);
        document.getElementById('cb-title').textContent = title;
    }

    getJetColor(v, min, max) {
        let t = (v - min) / (max - min || 1e-6);
        t = Math.min(Math.max(t, 0), 1.0);
        
        // Exact 2B.2 Heatmap Stops
        const stops = [
            { t: 0.0, r: 0.23, g: 0.51, b: 0.96 }, // #3b82f6 (Blue)
            { t: 0.2, r: 0.02, g: 0.71, b: 0.83 }, // #06b6d4 (Cyan)
            { t: 0.4, r: 0.06, g: 0.73, b: 0.51 }, // #10b981 (Green)
            { t: 0.6, r: 0.98, g: 0.80, b: 0.08 }, // #facc15 (Yellow)
            { t: 0.8, r: 0.98, g: 0.45, b: 0.09 }, // #f97316 (Orange)
            { t: 1.0, r: 0.94, g: 0.27, b: 0.27 }  // #ef4444 (Red)
        ];
        
        for (let i = 0; i < stops.length - 1; i++) {
            const s1 = stops[i];
            const s2 = stops[i+1];
            if (t >= s1.t && t <= s2.t) {
                const factor = (t - s1.t) / (s2.t - s1.t);
                return {
                    r: s1.r + factor * (s2.r - s1.r),
                    g: s1.g + factor * (s2.g - s1.g),
                    b: s1.b + factor * (s2.b - s1.b)
                };
            }
        }
        return stops[stops.length-1];
    }

    evaluateDisplacement(patch, u, v, disp) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        
        const spanU = this.nurbs.findSpan(nU - 1, p, u, U);
        const spanV = this.nurbs.findSpan(nV - 1, q, v, V);
        const dersU = this.nurbs.basisFuns(spanU, u, p, U);
        const dersV = this.nurbs.basisFuns(spanV, v, q, V);

        let ux = 0, uy = 0, W = 0;
        for (let i = 0; i <= p; i++) {
            for (let j = 0; j <= q; j++) {
                const cpI = spanU - p + i;
                const cpJ = spanV - q + j;
                const val = dersU[i] * dersV[j] * weights[cpI][cpJ];
                const idx = (cpI * nV + cpJ) * 2;
                ux += val * disp[idx];
                uy += val * disp[idx + 1];
                W += val;
            }
        }
        return { x: ux / W, y: uy / W };
    }

    renderControlNet(patch) {
        if (this.cpNet) this.scene.remove(this.cpNet);
        const lineMaterial = new THREE.LineBasicMaterial({ color: 0xf59e0b, transparent: true, opacity: 0.4 });
        const sphereMaterial = new THREE.MeshPhongMaterial({ color: 0xf59e0b, shininess: 100 });
        const sphereGeo = new THREE.SphereGeometry(0.04, 8, 8);
        
        const group = new THREE.Group();

        const cp = patch.controlPoints;
        const nU = cp.length;
        const nV = cp[0].length;

        // 1. Render Net Lines
        // U lines
        for (let i = 0; i < nU; i++) {
            const points = [];
            for (let j = 0; j < nV; j++) points.push(new THREE.Vector3(cp[i][j].x, cp[i][j].y, cp[i][j].z));
            const geo = new THREE.BufferGeometry().setFromPoints(points);
            group.add(new THREE.Line(geo, lineMaterial));
        }
        // V lines
        for (let j = 0; j < nV; j++) {
            const points = [];
            for (let i = 0; i < nU; i++) points.push(new THREE.Vector3(cp[i][j].x, cp[i][j].y, cp[i][j].z));
            const geo = new THREE.BufferGeometry().setFromPoints(points);
            group.add(new THREE.Line(geo, lineMaterial));
        }

        // 2. Render Spheres at Control Points
        for (let i = 0; i < nU; i++) {
            for (let j = 0; j < nV; j++) {
                const sphere = new THREE.Mesh(sphereGeo, sphereMaterial);
                sphere.position.set(cp[i][j].x, cp[i][j].y, cp[i][j].z);
                group.add(sphere);
            }
        }

        this.cpNet = group;
        this.cpNet.visible = this.showControlNet;
        this.scene.add(this.cpNet);

        this.renderBCs(patch);
    }

    renderBCs(patch) {
        if (this.bcGroup) this.scene.remove(this.bcGroup);
        this.bcGroup = new THREE.Group();
        
        const cp = patch.controlPoints;
        const nU = cp.length;
        const nV = cp[0].length;

        // 1. Symmetry BCs (Rollers)
        const rollerMaterial = new THREE.MeshPhongMaterial({ color: 0x94a3b8, transparent: true, opacity: 0.8 });
        const rollerGeo = new THREE.CylinderGeometry(0.05, 0.05, 0.02, 16);
        rollerGeo.rotateX(Math.PI / 2);

        // Bottom Edge (i=0) - Constrained in Y
        for (let j = 0; j < nV; j++) {
            const roller = new THREE.Mesh(rollerGeo, rollerMaterial);
            roller.position.set(cp[0][j].x, cp[0][j].y - 0.1, cp[0][j].z);
            this.bcGroup.add(roller);
        }

        // Left Edge (i=nU-1) - Constrained in X
        for (let j = 0; j < nV; j++) {
            const roller = new THREE.Mesh(rollerGeo, rollerMaterial);
            roller.position.set(cp[nU-1][j].x - 0.1, cp[nU-1][j].y, cp[nU-1][j].z);
            roller.rotateZ(Math.PI/2);
            this.bcGroup.add(roller);
        }

        // 2. Loads (Tension Arrows on Right Edge)
        const arrowDir = new THREE.Vector3(1, 0, 0);
        const arrowColor = 0xef4444; 
        const arrowLength = 0.4;

        const L = patch.params.L1;
        for (let i = 0; i < nU; i++) {
            const cp_node = cp[i][nV-1];
            if (cp_node.x > L - 0.1) {
                const origin = new THREE.Vector3(cp_node.x, cp_node.y, cp_node.z);
                const arrow = new THREE.ArrowHelper(arrowDir, origin, arrowLength, arrowColor, 0.15, 0.1);
                this.bcGroup.add(arrow);
            }
        }

        this.bcGroup.visible = this.showBCs !== false;
        this.scene.add(this.bcGroup);
    }

    setControlNetVisibility(visible) {
        this.showControlNet = visible;
        if (this.cpNet) this.cpNet.visible = visible;
    }

    setBCVisibility(visible) {
        this.showBCs = visible;
        if (this.bcGroup) this.bcGroup.visible = visible;
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        if (this.controls) this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
    renderECSWActiveMesh(patch, ecswData) {
        if (this.activeMeshGroup) this.scene.remove(this.activeMeshGroup);
        if (!ecswData || !ecswData.indices || !patch.elements) return;

        this.activeMeshGroup = new THREE.Group();
        const lineMaterial = new THREE.LineBasicMaterial({ 
            color: 0xfbbf24, 
            transparent: true, 
            opacity: 0.8,
            linewidth: 2
        });

        const points = [];
        const res = 2; // Low res for active mesh wireframe

        ecswData.indices.forEach(idx => {
            const element = patch.elements[idx];
            
            // Create 4 corners of the element for a wireframe box
            const corners = [
                {u: 0, v: 0}, {u: 1, v: 0},
                {u: 1, v: 0}, {u: 1, v: 1},
                {u: 1, v: 1}, {u: 0, v: 1},
                {u: 0, v: 1}, {u: 0, v: 0}
            ];

            corners.forEach(c => {
                const u = element.uRange[0] + c.u * (element.uRange[1] - element.uRange[0]);
                const v = element.vRange[0] + c.v * (element.vRange[1] - element.vRange[0]);
                const pos = this.nurbs.evaluateSurface(patch, u, v);
                points.push(new THREE.Vector3(pos.x, pos.y, pos.z + 0.012));
            });
        });

        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const segments = new THREE.LineSegments(geometry, lineMaterial);
        
        this.activeMeshGroup.add(segments);
        this.activeMeshGroup.visible = !!(window.app && window.app.params && window.app.params.showActiveMesh);
        this.scene.add(this.activeMeshGroup);
    }
}

window.VisualsEngine = VisualsEngine;
