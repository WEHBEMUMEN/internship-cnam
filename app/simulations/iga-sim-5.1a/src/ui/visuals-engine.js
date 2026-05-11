/**
 * Phase 5.1a - Visuals Engine (Three.js)
 */

class VisualsEngine {
    constructor(containerId, nurbs, solver) {
        this.container = document.getElementById(containerId);
        this.nurbs = nurbs;
        this.solver = solver;
        this.defScale = 50.0;
        this.showControlNet = true;
        this.viewMode = 'displacement'; // 'displacement' or 'stress'
        
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0f172a);
        
        this.camera = new THREE.PerspectiveCamera(45, this.container.clientWidth / this.container.clientHeight, 0.1, 1000);
        this.camera.position.set(2, 2, 10);
        this.camera.lookAt(2, 2, 0);

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
        this.controls.target.set(2, 2, 0);
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
                        if (s.vonMises > maxVal) maxVal = s.vonMises;
                    } catch (err) {
                        if (i === 0) console.error("[Visuals] Stress Eval Error:", err);
                    }
                } else if (this.viewMode === 'displacement') {
                    // Magnitude of displacement
                    const mag = Math.sqrt(u_disp.x**2 + u_disp.y**2);
                    vals[i] = mag;
                    if (mag > maxVal) maxVal = mag;
                }
            }
            positions.setXYZ(i, pos.x, pos.y, pos.z);
        }

        // 2. Map Colors (For Stress or Displacement Magnitude)
        if (maxVal > 0) {
            for (let i = 0; i < positions.count; i++) {
                const color = this.getJetColor(vals[i], 0, maxVal);
                colors[i * 3] = color.r;
                colors[i * 3 + 1] = color.g;
                colors[i * 3 + 2] = color.b;
            }
            geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        }

        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color: maxVal > 0 ? 0xffffff : 0x8b5cf6,
            vertexColors: maxVal > 0,
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

        this.renderControlNet(patch);
        
        // Update Color Bar
        const title = this.viewMode === 'stress' ? 'Von Mises (MPa)' : 'Displacement (mm)';
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
        let dv = max - min;
        if (dv === 0) dv = 1;
        const c = { r: 0, g: 0, b: 0 };
        const x = (v - min) / dv;
        if (x < 0.25) { c.r = 0; c.g = 4 * x; c.b = 1; }
        else if (x < 0.5) { c.r = 0; c.g = 1; c.b = 1 + 4 * (0.25 - x); }
        else if (x < 0.75) { c.r = 4 * (x - 0.5); c.g = 1; c.b = 0; }
        else { c.r = 1; c.g = 1 + 4 * (0.75 - x); c.b = 0; }
        return c;
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
        const rollerMaterial = new THREE.MeshPhongMaterial({ color: 0x94a3b8, transparent: true, opacity: 0.6 });
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

        // 2. Loads (Tension Arrows on Right Outer Edge j=nV-1, i <= nU/2)
        const arrowDir = new THREE.Vector3(1, 0, 0);
        const arrowColor = 0xef4444; // Red for force
        const nU_half = Math.floor(nU / 2);
        for (let i = 0; i <= nU_half; i++) {
            const origin = new THREE.Vector3(cp[i][nV-1].x, cp[i][nV-1].y, cp[i][nV-1].z);
            const arrow = new THREE.ArrowHelper(arrowDir, origin, 0.4, arrowColor, 0.15, 0.1);
            this.bcGroup.add(arrow);
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
}

window.VisualsEngine = VisualsEngine;
