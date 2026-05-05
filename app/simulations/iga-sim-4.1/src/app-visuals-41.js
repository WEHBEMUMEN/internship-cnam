/**
 * Phase 4.0 Visuals — Mesh Rendering & 3D Scene
 */

class TransientVisuals {
    constructor(app) {
        this.app = app;
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        
        this.deformedMesh = null;
        this.ghostMesh = null;
        this.wireframe = null;
        
        this.init();
    }

    init() {
        const container = document.getElementById('canvas-container');
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        container.appendChild(this.renderer.domElement);

        this.camera.position.set(2, 2, 5);
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        
        const ambient = new THREE.AmbientLight(0xffffff, 0.6);
        const direct = new THREE.DirectionalLight(0xffffff, 0.8);
        direct.position.set(5, 10, 5);
        this.scene.add(ambient, direct);

        window.addEventListener('resize', () => {
            this.camera.aspect = container.clientWidth / container.clientHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(container.clientWidth, container.clientHeight);
        });

        this.animate();
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    updateMesh(uDisp) {
        const res = 32;
        const posDef = [], colors = [];
        let maxF = 0;

        // 1. Calculate positions and find max displacement for scaling
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i/res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j/res, 0.9999);
                const st = this.app.fom.getSurfaceState(this.app.patch, u, v);
                let px = st.position.x, py = st.position.y;
                
                const d = this._interpDisp(u, v, st.denominator, uDisp);
                px += d.x; py += d.y;
                
                const mag = Math.sqrt(d.x*d.x + d.y*d.y);
                maxF = Math.max(maxF, mag);
                posDef.push(px, py, 0);
            }
        }

        // 2. Generate colors using Jet colormap
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i/res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j/res, 0.9999);
                const st = this.app.fom.getSurfaceState(this.app.patch, u, v);
                const d = this._interpDisp(u, v, st.denominator, uDisp);
                const f = maxF > 0 ? Math.sqrt(d.x*d.x + d.y*d.y) / maxF : 0;
                const [r, g, b] = this._jet(f);
                colors.push(r, g, b);
            }
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
            pa.copyArray(posDef);
            ca.copyArray(colors);
            pa.needsUpdate = true; ca.needsUpdate = true;
            this.deformedMesh.geometry.computeVertexNormals();
            this.wireframe.geometry.dispose();
            this.wireframe.geometry = new THREE.WireframeGeometry(this.deformedMesh.geometry);
        }
    }

    _interpDisp(u, v, denom, uDisp) {
        const nU = this.app.patch.controlPoints.length, nV = this.app.patch.controlPoints[0].length;
        const p = this.app.patch.p, q = this.app.patch.q;
        let dx = 0, dy = 0;
        
        // Use NURBS basis functions
        const spanU = window.nurbsUtils.findSpan(nU - 1, p, u, this.app.patch.U);
        const spanV = window.nurbsUtils.findSpan(nV - 1, q, v, this.app.patch.V);
        const Nu = window.nurbsUtils.dersBasisFuns(spanU, p, u, this.app.patch.U, 0)[0];
        const Nv = window.nurbsUtils.dersBasisFuns(spanV, q, v, this.app.patch.V, 0)[0];
        
        for (let i = 0; i <= p; i++) {
            for (let j = 0; j <= q; j++) {
                const idxI = spanU - p + i;
                const idxJ = spanV - q + j;
                const w = this.app.patch.weights[idxI][idxJ] || 1;
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
