/**
 * Phase 5.3 - Visuals Engine (Three.js)
 * High-performance side-by-side 2D/3D NURBS viewports rendering physical and reference domains.
 */

class VisualsEngine {
    constructor(containerId, isMapped = false) {
        this.container = document.getElementById(containerId);
        if (!this.container) {
            console.error(`[Visuals] Container #${containerId} not found!`);
            return;
        }

        this.isMapped = isMapped;
        this.defScale = 1.5;            // Scaled down displacement to look elegant
        this.showControlNet = false;
        this.showBCs = true;
        this.viewMode = 'stress';      // 'stress' or 'displacement'
        this.warpToPhysical = true;     // For mapped view: show morphed physical space or flat reference space
        
        // Scene setup
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0a0f1d); // Deep glassmorphic background
        
        // Camera setup
        this.camera = new THREE.PerspectiveCamera(45, this.container.clientWidth / this.container.clientHeight, 0.1, 1000);
        this.camera.position.set(5.0, 0.0, 14.0); // Perfect zoom for 10x2 cantilever beam
        this.camera.lookAt(5.0, 0.0, 0.0);

        // Renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.container.appendChild(this.renderer.domElement);

        // Controls
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false; // Keep it pure 2D for high clarity
        this.controls.screenSpacePanning = true;
        this.controls.mouseButtons = {
            LEFT: THREE.MOUSE.PAN,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.PAN
        };
        this.controls.target.set(5.0, 0.0, 0.0);
        this.controls.update();

        // Lighting
        const ambient = new THREE.AmbientLight(0xffffff, 0.9);
        this.scene.add(ambient);
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
        dirLight.position.set(5, 5, 10);
        this.scene.add(dirLight);

        // Object references
        this.mesh = null;
        this.wireframe = null;
        this.cpNetGroup = null;
        this.bcGroup = null;
        this.basisCache = null;
        this.cacheRes = 0;

        window.addEventListener('resize', () => this.onResize());
        this.animate();
    }

    onResize() {
        if (!this.container) return;
        this.camera.aspect = this.container.clientWidth / this.container.clientHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    /**
     * Cache NURBS shape functions for lightning-fast rendering
     */
    rebuildBasisCache(patch, resU, resV) {
        const count = (resU + 1) * (resV + 1);
        this.basisCache = new Array(count);
        this.cacheResU = resU;
        this.cacheResV = resV;

        const { p, q, U, V } = patch;
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const engine = new window.NURBS2D();

        let idx = 0;
        for (let j = 0; j <= resV; j++) {
            const v = j / resV;
            for (let i = 0; i <= resU; i++) {
                const u = i / resU;

                const spanU = engine.findSpan(nU - 1, p, u, U);
                const spanV = engine.findSpan(nV - 1, q, v, V);
                const dersU = engine.basisFuns(spanU, u, p, U);
                const dersV = engine.basisFuns(spanV, v, q, V);

                const contributions = [];
                for (let iu = 0; iu <= p; iu++) {
                    for (let jv = 0; jv <= q; jv++) {
                        const cpI = spanU - p + iu;
                        const cpJ = spanV - q + jv;
                        contributions.push({
                            cpI, cpJ,
                            basis: dersU[iu] * dersV[jv]
                        });
                    }
                }
                this.basisCache[idx++] = { u, v, contributions };
            }
        }
    }

    /**
     * Map physical coordinates via Pullback Mapping Φ(x_hat; mu, r)
     */
    evaluatePullbackMapping(hx, hy, mu, r, H) {
        const dist = Math.abs(hx - mu);
        const sigma = Math.max(r, 0.4);
        const disp = r * Math.exp(-0.5 * Math.pow(dist / sigma, 2));

        // Morphed coordinates: bottom notch curves up, top notch curves down
        let y = hy;
        const scaleFactor = 1.0 - (2.0 * disp / H);
        y = hy * scaleFactor;

        return { x: hx, y };
    }

    /**
     * Update the viewport with the solved state
     */
    update(solver, patch, displacement = null) {
        if (this.mesh) this.scene.remove(this.mesh);
        if (this.wireframe) this.scene.remove(this.wireframe);
        if (this.cpNetGroup) this.scene.remove(this.cpNetGroup);
        if (this.bcGroup) this.scene.remove(this.bcGroup);

        const resU = 40;
        const resV = 16;
        
        if (!this.basisCache || this.cacheResU !== resU || this.cacheResV !== resV) {
            this.rebuildBasisCache(patch, resU, resV);
        }

        const geometry = new THREE.PlaneGeometry(1, 1, resU, resV);
        const positions = geometry.attributes.position;
        const colors = new Float32Array(positions.count * 3);

        const cp = patch.controlPoints;
        const weights = patch.weights;
        const nV = cp[0].length;
        const nNodes = cp.length * cp[0].length;

        // Evaluate max values for colormap
        let maxVal = 0.0;
        const vals = new Float32Array(positions.count);

        if (displacement) {
            if (this.viewMode === 'displacement') {
                for (let i = 0; i < displacement.length; i += 2) {
                    const mag = Math.sqrt(displacement[i]**2 + displacement[i+1]**2);
                    if (mag > maxVal) maxVal = mag;
                }
            } else {
                // Von Mises stress sampling
                // Sample 5 points to scale colormap
                let tempMax = 0;
                const stressSolver = new window.IGA2DSolver(solver.engine);
                for (let i = 0; i < positions.count; i += 8) {
                    const cache = this.basisCache[i];
                    try {
                        const s = stressSolver.getNumericalStress(patch, displacement, cache.u, cache.v, solver.E, solver.nu);
                        if (s.vonMises > tempMax) tempMax = s.vonMises;
                    } catch(e) {}
                }
                maxVal = tempMax || 1.0;
            }
        }
        if (maxVal === 0) maxVal = 1.0;

        const stressSolver = new window.IGA2DSolver(solver.engine);

        // Build deformed coordinates and contours
        for (let idx = 0; idx < positions.count; idx++) {
            const cache = this.basisCache[idx];
            let x = 0, y = 0, z = 0, W = 0;
            let ux = 0, uy = 0;

            for (let i = 0; i < cache.contributions.length; i++) {
                const c = cache.contributions[i];
                const weight = weights[c.cpI][c.cpJ];
                const val = c.basis * weight;
                
                x += val * cp[c.cpI][c.cpJ].x;
                y += val * cp[c.cpI][c.cpJ].y;
                z += val * cp[c.cpI][c.cpJ].z;
                
                if (displacement) {
                    const dIdx = (c.cpI * nV + c.cpJ) * 2;
                    ux += val * displacement[dIdx];
                    uy += val * displacement[dIdx + 1];
                }
                W += val;
            }

            let posX = x / W;
            let posY = y / W;

            // In Mapped viewport, warp coordinates forward to physical space if requested
            if (this.isMapped) {
                if (this.warpToPhysical) {
                    const morphed = this.evaluatePullbackMapping(posX, posY, solver.mu, solver.r, solver.H);
                    posX = morphed.x;
                    posY = morphed.y;
                } else {
                    // Force strictly flat reference configuration grid
                    posX = (idx % (resU + 1)) / resU * solver.L;
                    posY = ((Math.floor(idx / (resU + 1))) / resV - 0.5) * solver.H;
                }
            }

            // Apply displacement deformation scale
            const defX = posX + (ux / W) * this.defScale;
            const defY = posY + (uy / W) * this.defScale;
            
            positions.setXYZ(idx, defX, defY, z / W);

            if (displacement) {
                if (this.viewMode === 'stress') {
                    try {
                        const s = stressSolver.getNumericalStress(patch, displacement, cache.u, cache.v, solver.E, solver.nu);
                        vals[idx] = s.vonMises;
                    } catch(e) {
                        vals[idx] = 0;
                    }
                } else {
                    vals[idx] = Math.sqrt(ux*ux + uy*uy);
                }
            }
        }

        // Apply dynamic colormap
        for (let i = 0; i < positions.count; i++) {
            const color = this.getContourColor(vals[i], 0, maxVal);
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
        }

        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        geometry.computeVertexNormals();

        // Sleek Material with beautiful specular highlights
        const material = new THREE.MeshPhongMaterial({
            color: 0xffffff,
            vertexColors: true,
            side: THREE.DoubleSide,
            shininess: 80,
            specular: 0x444444,
            transparent: true,
            opacity: 0.95
        });

        this.mesh = new THREE.Mesh(geometry, material);
        this.scene.add(this.mesh);

        // Thin elegant gridlines
        const wireMaterial = new THREE.MeshBasicMaterial({
            color: this.isMapped ? 0x8b5cf6 : 0x06b6d4, // Violet for Mapped, Cyan for Naive
            wireframe: true,
            transparent: true,
            opacity: 0.25
        });
        this.wireframe = new THREE.Mesh(geometry, wireMaterial);
        this.scene.add(this.wireframe);

        // Optional rendering overlays
        if (this.showControlNet) {
            this.renderControlNet(patch, displacement, solver);
        }
        if (this.showBCs) {
            this.renderBCs(solver, patch);
        }

        this.updateLegend(0, maxVal);
    }

    /**
     * Render the B-Spline control mesh overlay
     */
    renderControlNet(patch, displacement, solver) {
        this.cpNetGroup = new THREE.Group();
        const cp = patch.controlPoints;
        const nU = cp.length;
        const nV = cp[0].length;

        const sphereGeo = new THREE.SphereGeometry(0.08, 16, 16);
        const sphereMat = new THREE.MeshPhongMaterial({
            color: this.isMapped ? 0xc084fc : 0x22d3ee,
            emissive: this.isMapped ? 0x3b0764 : 0x083344,
            shininess: 100
        });

        const points = [];

        for (let i = 0; i < nU; i++) {
            points[i] = [];
            for (let j = 0; j < nV; j++) {
                let x = cp[i][j].x;
                let y = cp[i][j].y;

                if (this.isMapped && !this.warpToPhysical) {
                    x = (i / (nU - 1)) * solver.L;
                    y = (j / (nV - 1) - 0.5) * solver.H;
                }

                if (displacement) {
                    const dIdx = (i * nV + j) * 2;
                    x += displacement[dIdx] * this.defScale;
                    y += displacement[dIdx + 1] * this.defScale;
                }

                const posVec = new THREE.Vector3(x, y, cp[i][j].z);
                points[i][j] = posVec;

                // Sphere at control point
                const sphere = new THREE.Mesh(sphereGeo, sphereMat);
                sphere.position.copy(posVec);
                this.cpNetGroup.add(sphere);
            }
        }

        // Draw grid lines connecting control points
        const lineMat = new THREE.LineBasicMaterial({
            color: this.isMapped ? 0xa855f7 : 0x06b6d4,
            transparent: true,
            opacity: 0.5
        });

        // U-direction lines
        for (let j = 0; j < nV; j++) {
            const linePoints = [];
            for (let i = 0; i < nU; i++) {
                linePoints.push(points[i][j]);
            }
            const lineGeo = new THREE.BufferGeometry().setFromPoints(linePoints);
            const line = new THREE.Line(lineGeo, lineMat);
            this.cpNetGroup.add(line);
        }

        // V-direction lines
        for (let i = 0; i < nU; i++) {
            const linePoints = [];
            for (let j = 0; j < nV; j++) {
                linePoints.push(points[i][j]);
            }
            const lineGeo = new THREE.BufferGeometry().setFromPoints(linePoints);
            const line = new THREE.Line(lineGeo, lineMat);
            this.cpNetGroup.add(line);
        }

        this.scene.add(this.cpNetGroup);
    }

    /**
     * Render Dirichlet support boundary locks at the left clamped edge
     */
    renderBCs(solver, patch) {
        this.bcGroup = new THREE.Group();
        const cp = patch.controlPoints;
        const nV = cp[0].length;

        const triangleGeo = new THREE.ConeGeometry(0.15, 0.3, 4);
        const triangleMat = new THREE.MeshBasicMaterial({ color: 0xef4444 });

        for (let j = 0; j < nV; j++) {
            let x = cp[0][j].x;
            let y = cp[0][j].y;

            if (this.isMapped && !this.warpToPhysical) {
                x = 0.0;
                y = (j / (nV - 1) - 0.5) * solver.H;
            }

            const anchor = new THREE.Mesh(triangleGeo, triangleMat);
            anchor.position.set(x - 0.25, y, 0);
            anchor.rotation.z = -Math.PI / 2; // Pointing to clamp the boundary
            this.bcGroup.add(anchor);
        }

        this.scene.add(this.bcGroup);
    }

    /**
     * Map value to a vibrant HSL-derived JET-like colormap
     */
    getContourColor(value, min, max) {
        let pct = (value - min) / (max - min);
        pct = Math.max(0, Math.min(1, pct)); // Clamp to 0..1

        // Smooth Viridis/Jet hybrid palette for modern aesthetics
        // Cyan (0.0) -> Purple/Violet (0.5) -> Intense Red (1.0)
        let r = 0, g = 0, b = 0;
        
        if (pct < 0.33) {
            // Cyan to Blue
            const t = pct / 0.33;
            r = 0;
            g = Math.round(180 * t);
            b = Math.round(180 + 75 * t);
        } else if (pct < 0.66) {
            // Blue to Violet
            const t = (pct - 0.33) / 0.33;
            r = Math.round(139 * t);
            g = Math.round(180 - 100 * t);
            b = Math.round(255 - 50 * t);
        } else {
            // Violet to Neon Orange/Red
            const t = (pct - 0.66) / 0.34;
            r = Math.round(139 + 116 * t);
            g = Math.round(80 - 60 * t);
            b = Math.round(205 - 180 * t);
        }

        return { r: r / 255, g: g / 255, b: b / 255 };
    }

    /**
     * Update the legend displays
     */
    updateLegend(min, max) {
        const prefix = this.isMapped ? 'mapped' : 'naive';
        const labelMin = document.getElementById(`${prefix}-legend-min`);
        const labelMax = document.getElementById(`${prefix}-legend-max`);
        const labelTitle = document.getElementById(`${prefix}-legend-title`);

        if (labelMin) labelMin.innerText = min.toFixed(1);
        if (labelMax) labelMax.innerText = max.toFixed(1);
        if (labelTitle) {
            labelTitle.innerText = this.viewMode === 'stress' ? 'von Mises Stress (MPa)' : 'Displacement Magnitude (mm)';
        }
    }
}

window.VisualsEngine = VisualsEngine;
