
// Verification Lab 2.3 | h-Refinement (Knot Insertion)
// Phase 2.3 | Computational Core

const engine = new NURBS2D();
let patch = NURBSPresets.generateSheet();

// --- Three.js Setup ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 2000);
camera.position.set(40, 40, 60);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const orbitControls = new THREE.OrbitControls(camera, renderer.domElement);
const transformControls = new THREE.TransformControls(camera, renderer.domElement);
scene.add(transformControls);

scene.add(new THREE.AmbientLight(0xffffff, 0.5));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(20, 100, 20);
scene.add(dLight);
scene.add(new THREE.GridHelper(100, 20, 0x334155, 0x1e293b));

// --- Objects ---
let surfaceMesh, wireframeOverlay;
let isoparmaGroup = new THREE.Group();
let gridGroup = new THREE.Group();
let pointMeshes = []; 
scene.add(gridGroup, isoparmaGroup);

// --- Core Logic ---

function syncUI() {
    document.getElementById('state-degree').textContent = `Degree: p=${patch.p}, q=${patch.q}`;
    document.getElementById('state-points').textContent = `Points: ${patch.controlPoints.length} x ${patch.controlPoints[0].length}`;
}

function fullRebuild() {
    transformControls.detach();
    pointMeshes = [];
    gridGroup.clear();
    createSurface();
    updateGrid();
    syncUI();
}

function createSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);
    isoparmaGroup.clear();
    
    const res = parseInt(document.getElementById('resolution-slider').value);
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            const p = state.position;
            positions.push(p.x, p.y, p.z);
            
            const detJ = engine.getJacobianDeterminant(patch, u, v);
            const val = Math.min(Math.max(detJ / 100, 0.0), 1.0);
            colors.push(0.1, 0.5, 0.2 + 0.8 * val); 
        }
    }

    const indices = [];
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            const a = i * (res + 1) + j;
            const b = (i + 1) * (res + 1) + j;
            const c = (i + 1) * (res + 1) + (j + 1);
            const d = i * (res + 1) + (j + 1);
            indices.push(a, b, d, b, c, d);
        }
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshStandardMaterial({
        vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.8
    });
    
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);

    if (document.getElementById('show-wireframe').checked) {
        wireframeOverlay = new THREE.LineSegments(new THREE.WireframeGeometry(geometry), new THREE.LineBasicMaterial({ color: 0x34d399, transparent: true, opacity: 0.2 }));
        scene.add(wireframeOverlay);
    }
}

function updateGrid() {
    gridGroup.clear();
    const lineMat = new THREE.LineBasicMaterial({ color: 0x475569, transparent: true, opacity: 0.4 });
    const n = patch.controlPoints.length;
    const m = patch.controlPoints[0].length;

    if (pointMeshes.length === 0) {
        const sphereGeom = new THREE.SphereGeometry(0.4, 12, 12);
        const sphereMat = new THREE.MeshBasicMaterial({ color: 0x10b981 }); // Greenish for h-refine
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < m; j++) {
                const sp = new THREE.Mesh(sphereGeom, sphereMat.clone());
                sp.userData = { i, j };
                pointMeshes.push(sp);
            }
        }
    }

    if (document.getElementById('show-grid').checked) {
        pointMeshes.forEach(sp => {
            const { i, j } = sp.userData;
            const cp = patch.controlPoints[i][j];
            if (transformControls.object !== sp) sp.position.set(cp.x, cp.y, cp.z);
            gridGroup.add(sp);
        });

        // Lattice lines
        for (let i = 0; i < n; i++) {
            const pts = [];
            for (let j = 0; j < m; j++) {
                const cp = patch.controlPoints[i][j];
                pts.push(cp.x, cp.y, cp.z);
            }
            gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(pts, 3)), lineMat));
        }
        for (let j = 0; j < m; j++) {
            const pts = [];
            for (let i = 0; i < n; i++) {
                const cp = patch.controlPoints[i][j];
                pts.push(cp.x, cp.y, cp.z);
            }
            gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(pts, 3)), lineMat));
        }
    }
}

// Interaction
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

function onPointerDown(event) {
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(pointMeshes);
    if (intersects.length > 0) {
        const s = intersects[0].object;
        transformControls.attach(s);
        pointMeshes.forEach(p => p.material.color.set(0x10b981)); 
        s.material.color.set(0xffffff); 
    } else if (!transformControls.dragging) {
        transformControls.detach();
        pointMeshes.forEach(p => p.material.color.set(0x10b981));
    }
}

transformControls.addEventListener('dragging-changed', (e) => orbitControls.enabled = !e.value);
transformControls.addEventListener('objectChange', () => {
    const obj = transformControls.object;
    if (obj) {
        const { i, j } = obj.userData;
        patch.controlPoints[i][j] = { x: obj.position.x, y: obj.position.y, z: obj.position.z };
        createSurface();
        updateGrid();
    }
});

// Operations
document.getElementById('refine-h').onclick = () => {
    engine.insertKnotU(patch, 0.5);
    engine.insertKnotV(patch, 0.5);
    fullRebuild();
};

document.getElementById('reset-surface').onclick = () => {
    patch = NURBSPresets.generateSheet();
    fullRebuild();
};

document.getElementById('resolution-slider').oninput = createSurface;
document.getElementById('show-grid').onchange = updateGrid;
document.getElementById('show-wireframe').onchange = createSurface;
renderer.domElement.addEventListener('pointerdown', onPointerDown);

window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

function animate() {
    requestAnimationFrame(animate);
    orbitControls.update();
    renderer.render(scene, camera);
}

fullRebuild();
animate();
