/**
 * Phase 4.1 Visuals — Mesh Rendering & 3D Scene
 */

class TransientVisuals {
    constructor(app) {
        this.app = app;
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        
        this.deformedMesh = null;
        this.wireframe = null;
        
        this.init();
    }

    init() {
        const container = document.getElementById('canvas-container');
        if (!container) return;
        
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        container.appendChild(this.renderer.domElement);

        this.camera.position.set(10, 10, 30);
        this.camera.lookAt(0, 0, 0);
        
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        const ambient = new THREE.AmbientLight(0xffffff, 0.6);
        const direct = new THREE.DirectionalLight(0xffffff, 0.8);
        direct.position.set(5, 10, 5);
        this.scene.add(ambient, direct);

        window.addEventListener('resize', () => {
            if (!container) return;
            this.camera.aspect = container.clientWidth / container.clientHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(container.clientWidth, container.clientHeight);
        });

        this.animate();
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        if (this.controls) this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    updateMesh(uDisp) {
        const res = 24; 
        const posDef = [], colors = [];
        let maxDisp = 0;

        const engine = this.app.engine;
        const patch = this.app.patch;

        const states = [];
        for (let i = 0; i <= res; i++) {
            const xi = i / res;
            for (let j = 0; j <= res; j++) {
                const eta = j / res;
                const st = engine.getSurfaceState(patch, xi, eta);
                const d = this._interpDisp(xi, eta, st.denominator, uDisp);
                const mag = Math.sqrt(d.x*d.x + d.y*d.y);
                maxDisp = Math.max(maxDisp, mag);
                states.push({ st, d, mag });
            }
        }

        for (const s of states) {
            posDef.push(s.st.position.x + s.d.x, s.st.position.y + s.d.y, s.st.position.z);
            const f = maxDisp > 0 ? s.mag / maxDisp : 0;
            const [r, g, b] = this._jet(f);
            colors.push(r, g, b);
        }

        if (!this.deformedMesh) {
            const idx = [];
            for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
                const a = i*(res+1)+j, b = (i+1)*(res+1)+j, c = (i+1)*(res+1)+j+1, d = i*(res+1)+j+1;
                idx.push(a,b,d, b,c,d);
            }
            const g = new THREE.BufferGeometry();
            g.setIndex(idx);
            g.setAttribute('position', new THREE.Float32BufferAttribute(posDef, 3));
            g.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
            
            this.deformedMesh = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: 0.5, metalness: 0.2 }));
            this.wireframe = new THREE.LineSegments(new THREE.WireframeGeometry(g), new THREE.LineBasicMaterial({ color: 0x0ea5e9, transparent: true, opacity: 0.2 }));
            
            this.scene.add(this.deformedMesh, this.wireframe);
        } else {
            const pa = this.deformedMesh.geometry.getAttribute('position');
            const ca = this.deformedMesh.geometry.getAttribute('color');
            pa.copyArray(new Float32Array(posDef));
            ca.copyArray(new Float32Array(colors));
            pa.needsUpdate = true; ca.needsUpdate = true;
            this.deformedMesh.geometry.computeVertexNormals();
            this.wireframe.geometry.dispose();
            this.wireframe.geometry = new THREE.WireframeGeometry(this.deformedMesh.geometry);
        }
    }

    _interpDisp(xi, eta, denom, uDisp) {
        const patch = this.app.patch;
        const engine = this.app.engine;
        const { p, q, U, V, weights, controlPoints } = patch;
        const nV = controlPoints[0].length;
        const nU = controlPoints.length;

        const spanU = engine.findSpan(nU - 1, p, xi, U);
        const spanV = engine.findSpan(nV - 1, q, eta, V);
        const Nu = engine.basisFuns(spanU, xi, p, U);
        const Nv = engine.basisFuns(spanV, eta, q, V);

        let dx = 0, dy = 0;
        for (let i = 0; i <= p; i++) {
            for (let j = 0; j <= q; j++) {
                const idxI = spanU - p + i;
                const idxJ = spanV - q + j;
                const w = weights[idxI][idxJ];
                const R = (Nu[i] * Nv[j] * w) / denom;
                const dofIdx = (idxI * nV + idxJ) * 2;
                dx += R * uDisp[dofIdx];
                dy += R * uDisp[dofIdx + 1];
            }
        }
        return { x: dx, y: dy };
    }

    _jet(v) {
        const r = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * v - 3)));
        const g = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * v - 2)));
        const b = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * v - 1)));
        return [r, g, b];
    }
}
