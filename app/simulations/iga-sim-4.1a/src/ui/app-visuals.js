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
        this.markersGroup = new THREE.Group();
        this.forceArrow = null;
        this.scene.add(this.markersGroup);
        
        this.defScale = 10;
        this.showCP = false;
        this.showActiveElements = false;
        this.controlPointsGroup = new THREE.Group();
        this.activeElementsGroup = new THREE.Group();
        this.scene.add(this.controlPointsGroup, this.activeElementsGroup);
        this.init();
    }

    init() {
        const container = document.getElementById('canvas-container');
        if (!container) return;
        
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        container.appendChild(this.renderer.domElement);

        this.camera.position.set(5, 1, 25);
        this.camera.lookAt(5, 1, 0);
        
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.target.set(5, 1, 0);
        this.controls.enableRotate = true; 
        this.controls.mouseButtons = {
            LEFT: THREE.MOUSE.PAN,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.ROTATE
        };
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

    drawMarkers() {
        this.markersGroup.clear();
        const dyn = this.app.dyn;
        if (!dyn || !dyn.patch) return;

        const patch = dyn.patch;
        const nV = patch.controlPoints[0].length;
        const nU = patch.controlPoints.length;

        const coneGeo = new THREE.ConeGeometry(0.15, 0.35, 4); 
        const coneMat = new THREE.MeshBasicMaterial({ color: 0x3b82f6, depthTest: false });

        const constraint = patch.constraints?.[0] || { side: 'left' };
        const side = constraint.side;

        for (let i = 0; i < nU; i++) {
            for (let j = 0; j < nV; j++) {
                
                let isFixed = false;
                if (side === 'left' && i === 0) isFixed = true;
                else if (side === 'right' && i === nU - 1) isFixed = true;
                else if (side === 'bottom' && j === 0) isFixed = true;
                else if (side === 'top' && j === nV - 1) isFixed = true;

                if (isFixed) {
                    const cp = patch.controlPoints[i][j];
                    const cone = new THREE.Mesh(coneGeo, coneMat);
                    cone.renderOrder = 999;
                    
                    let offsetX = 0, offsetY = 0, rotZ = 0;
                    if (side === 'left') { offsetX = -0.3; rotZ = -Math.PI/2; }
                    else if (side === 'right') { offsetX = 0.3; rotZ = Math.PI/2; }
                    else if (side === 'bottom') { offsetY = -0.3; rotZ = 0; }
                    else if (side === 'top') { offsetY = 0.3; rotZ = Math.PI; }

                    cone.position.set(cp.x + offsetX, cp.y + offsetY, cp.z);
                    cone.rotation.z = rotZ;
                    this.markersGroup.add(cone);
                }
            }
        }
    }

    drawControlPoints(uDisp) {
        this.controlPointsGroup.clear();
        if (!this.showCP || !uDisp) return;

        const patch = this.app.patch;
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const scale = this.defScale;

        const sphereGeo = new THREE.SphereGeometry(0.1, 8, 8);
        const sphereMat = new THREE.MeshBasicMaterial({ color: 0xffffff, depthTest: false });

        for (let i = 0; i < nU; i++) {
            for (let j = 0; j < nV; j++) {
                const cp = patch.controlPoints[i][j];
                const dofIdx = (i * nV + j) * 2;
                const sphere = new THREE.Mesh(sphereGeo, sphereMat);
                sphere.position.set(
                    cp.x + uDisp[dofIdx] * scale,
                    cp.y + uDisp[dofIdx+1] * scale,
                    cp.z
                );
                sphere.renderOrder = 2000;
                this.controlPointsGroup.add(sphere);
            }
        }
    }

    updateActiveElements() {
        this.activeElementsGroup.clear();
        if (!this.showActiveElements || !this.app.trainer || !this.app.trainer.package.ecsw) return;

        const patch = this.app.patch;
        const { U, V } = patch;
        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const nV = uniqueV.length - 1;

        const mat = new THREE.LineBasicMaterial({ color: 0x10b981, linewidth: 2, depthTest: false });

        this.app.trainer.package.ecsw.indices.forEach(idx => {
            const i = Math.floor(idx / nV);
            const j = idx % nV;
            
            const uMin = uniqueU[i], uMax = uniqueU[i+1];
            const vMin = uniqueV[j], vMax = uniqueV[j+1];

            const pts = [];
            const res = 4;
            for(let k=0; k<=res; k++) pts.push(this.app.engine.getSurfaceState(patch, uMin + (uMax-uMin)*(k/res), vMin).position);
            for(let k=0; k<=res; k++) pts.push(this.app.engine.getSurfaceState(patch, uMax, vMin + (vMax-vMin)*(k/res)).position);
            for(let k=0; k<=res; k++) pts.push(this.app.engine.getSurfaceState(patch, uMax - (uMax-uMin)*(k/res), vMax).position);
            for(let k=0; k<=res; k++) pts.push(this.app.engine.getSurfaceState(patch, uMin, vMax - (vMax-vMin)*(k/res)).position);

            const geometry = new THREE.BufferGeometry().setFromPoints(pts.map(p => new THREE.Vector3(p.x, p.y, p.z + 0.1)));
            const line = new THREE.Line(geometry, mat);
            line.renderOrder = 3000;
            this.activeElementsGroup.add(line);
        });
    }

    updateForceMarkers(t) {
        if (!this.forceArrowsGroup) {
            this.forceArrowsGroup = new THREE.Group();
            this.scene.add(this.forceArrowsGroup);
        }
        this.forceArrowsGroup.clear();

        const load = this.app.dyn.load;
        if (!load) return;

        const timeFactor = load.timeFunction ? load.timeFunction(t) : 1.0;
        const engine = this.app.engine;
        const patch = this.app.patch;

        const renderArrow = (pos, magX, magY) => {
            const mag = Math.sqrt(magX*magX + magY*magY);
            if (mag < 1e-3) return;
            const dir = new THREE.Vector3(magX, magY, 0).normalize();
            const arrow = new THREE.ArrowHelper(dir, pos, mag * 0.5, 0xf43f5e, 0.4, 0.2);
            arrow.renderOrder = 1000;
            this.forceArrowsGroup.add(arrow);
        };

        if (load.type === 'distributed') {
            const res = 5; 
            for (let k = 0; k <= res; k++) {
                const eta = k / res;
                const st = engine.getSurfaceState(patch, 1.0, eta); 
                const pos = new THREE.Vector3(st.position.x, st.position.y, st.position.z);
                renderArrow(pos, load.magnitude.x * timeFactor, load.magnitude.y * timeFactor);
            }
        } else if (load.type === 'point') {
            const st = engine.getSurfaceState(patch, load.xi, load.eta);
            const pos = new THREE.Vector3(st.position.x, st.position.y, st.position.z);
            renderArrow(pos, load.magnitude.x * timeFactor, load.magnitude.y * timeFactor);
        }
    }

    updateMesh(uDisp, t = 0) {
        if (this.markersGroup.children.length === 0) this.drawMarkers();
        this.updateForceMarkers(t);
        this.drawControlPoints(uDisp);
        this.updateActiveElements();
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

        const scale = this.defScale;
        for (const s of states) {
            posDef.push(
                s.st.position.x + s.d.x * scale, 
                s.st.position.y + s.d.y * scale, 
                s.st.position.z
            );
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
