/**
 * Phase 5.1a - Visuals Engine (Three.js)
 */

class VisualsEngine {
    constructor(containerId, nurbs) {
        this.container = document.getElementById(containerId);
        this.nurbs = nurbs;
        this.defScale = 50.0;
        
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
            RIGHT: THREE.MOUSE.ROTATE
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

        const geometry = new THREE.PlaneGeometry(1, 1, 32, 32);
        const positions = geometry.attributes.position;
        
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        // Evaluate surface
        for (let i = 0; i < positions.count; i++) {
            const u = geometry.attributes.uv.getX(i);
            const v = geometry.attributes.uv.getY(i);
            
            // Basic NURBS evaluation (using injected engine)
            const pos = this.nurbs.evaluateSurface(patch, u, v);
            
            // Apply displacement if provided
            if (displacement) {
                const u_disp = this.evaluateDisplacement(patch, u, v, displacement);
                pos.x += u_disp.x * this.defScale;
                pos.y += u_disp.y * this.defScale;
            }

            positions.setXYZ(i, pos.x, pos.y, pos.z);
        }

        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color: 0x8b5cf6,
            side: THREE.DoubleSide,
            flatShading: false,
            transparent: true,
            opacity: 0.8,
            shininess: 100
        });

        this.mesh = new THREE.Mesh(geometry, material);
        this.scene.add(this.mesh);

        const wireMaterial = new THREE.MeshBasicMaterial({ color: 0x475569, wireframe: true, transparent: true, opacity: 0.3 });
        this.wireframe = new THREE.Mesh(geometry, wireMaterial);
        this.scene.add(this.wireframe);

        this.renderControlNet(patch);
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
        const material = new THREE.LineBasicMaterial({ color: 0xf59e0b, transparent: true, opacity: 0.4 });
        const group = new THREE.Group();

        const cp = patch.controlPoints;
        const nU = cp.length;
        const nV = cp[0].length;

        // U lines
        for (let i = 0; i < nU; i++) {
            const points = [];
            for (let j = 0; j < nV; j++) points.push(new THREE.Vector3(cp[i][j].x, cp[i][j].y, cp[i][j].z));
            const geo = new THREE.BufferGeometry().setFromPoints(points);
            group.add(new THREE.Line(geo, material));
        }
        // V lines
        for (let j = 0; j < nV; j++) {
            const points = [];
            for (let i = 0; i < nU; i++) points.push(new THREE.Vector3(cp[i][j].x, cp[i][j].y, cp[i][j].z));
            const geo = new THREE.BufferGeometry().setFromPoints(points);
            group.add(new THREE.Line(geo, material));
        }

        this.cpNet = group;
        this.scene.add(this.cpNet);
    }

    setControlNetVisibility(visible) {
        if (this.cpNet) this.cpNet.visible = visible;
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.renderer.render(this.scene, this.camera);
    }
}

window.VisualsEngine = VisualsEngine;
