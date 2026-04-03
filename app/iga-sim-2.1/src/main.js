
// Verification Lab 2.1 | Interactive Surface Mapping
// Phase 2.1 | Computational Core

const engine = new NURBS2D();

// --- Patch Setup (Quadratic 4x4) ---
const patch = {
    p: 2, q: 2,
    U: [0, 0, 0, 0.5, 1, 1, 1], // Open knot vector
    V: [0, 0, 0, 0.5, 1, 1, 1],
    weights: Array(4).fill().map(() => Array(4).fill(1)),
    controlPoints: []
};

function initDefaultPoints() {
    patch.controlPoints = [];
    for (let i = 0; i < 4; i++) {
        patch.controlPoints[i] = [];
        for (let j = 0; j < 4; j++) {
            patch.controlPoints[i][j] = {
                x: (i - 1.5) * 10,
                y: (j - 1.5) * 10,
                z: 0
            };
        }
    }
    // Default "Sheet" curvature
    patch.controlPoints[0][0].z = 5; patch.controlPoints[0][3].z = -5;
    patch.controlPoints[3][0].z = -5; patch.controlPoints[3][3].z = 5;
    patch.controlPoints[1][1].z = 15; patch.controlPoints[2][2].z = 15;
    patch.controlPoints[1][2].z = -10; patch.controlPoints[2][1].z = -10;
}
initDefaultPoints();

// --- Three.js Setup ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const orbitControls = new THREE.OrbitControls(camera, renderer.domElement);
const transformControls = new THREE.TransformControls(camera, renderer.domElement);
transformControls.setMode('translate');
transformControls.setSize(1.5); // Larger gizmo for better visibility
scene.add(transformControls);

const ambientLight = new THREE.AmbientLight(0xffffff, 0.3);
const directionalLight = new THREE.DirectionalLight(0xffffff, 1.0);
directionalLight.position.set(20, 40, 20);
scene.add(ambientLight, directionalLight);

const pointLight = new THREE.PointLight(0x3b82f6, 1, 100);
pointLight.position.set(0, 20, 10);
scene.add(pointLight);

camera.position.set(30,30,50);
orbitControls.update();

// --- Objects ---
let surfaceMesh, wireframeOverlay;
let isoparmaGroup = new THREE.Group();
let gridGroup = new THREE.Group();
let pointMeshes = []; // 1D array of 16 spheres
scene.add(gridGroup, isoparmaGroup);

// --- Core Logic ---

function createSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);
    isoparmaGroup.clear();
    
    const resolution = parseInt(document.getElementById('resolution-slider').value);
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    for (let i = 0; i <= resolution; i++) {
        const u = i / resolution;
        for (let j = 0; j <= resolution; j++) {
            const v = j / resolution;
            const state = engine.getSurfaceState(patch, u, v);
            const p = state.position;
            positions.push(p.x, p.y, p.z);
            
            const detJ = engine.getJacobianDeterminant(patch, u, v);
            const val = Math.min(detJ / 100, 1.0);
            colors.push(0.1, 0.4, 0.8 * val + 0.2);
        }
    }

    const indices = [];
    for (let i = 0; i < resolution; i++) {
        for (let j = 0; j < resolution; j++) {
            const a = i * (resolution + 1) + j;
            const b = (i + 1) * (resolution + 1) + j;
            const c = (i + 1) * (resolution + 1) + (j + 1);
            const d = i * (resolution + 1) + (j + 1);
            indices.push(a, b, d, b, c, d);
        }
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshStandardMaterial({
        vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.85, roughness: 0.3, metalness: 0.2
    });
    
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);

    if (document.getElementById('show-wireframe').checked) {
        wireframeOverlay = new THREE.LineSegments(new THREE.WireframeGeometry(geometry), new THREE.LineBasicMaterial({ color: 0x60a5fa, transparent: true, opacity: 0.3 }));
        scene.add(wireframeOverlay);
    }

    // Isoparms
    const isoCount = 10;
    const isoMat = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.4 });
    for (let k = 0; k <= isoCount; k++) {
        const t = k / isoCount;
        const ptsU = [], ptsV = [];
        for (let s = 0; s <= resolution; s++) {
            const step = s / resolution;
            const pU = engine.evaluateSurface(patch, t, step);
            const pV = engine.evaluateSurface(patch, step, t);
            ptsU.push(pU.x, pU.y, pU.z); ptsV.push(pV.x, pV.y, pV.z);
        }
        isoparmaGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(ptsU, 3)), isoMat));
        isoparmaGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(ptsV, 3)), isoMat));
    }
}

function updateGrid() {
    gridGroup.clear();
    if (!document.getElementById('show-grid').checked) return;

    const lineMat = new THREE.LineBasicMaterial({ color: 0x475569, transparent: true, opacity: 0.6 });

    // Sync sphere meshes OR create them
    if (pointMeshes.length === 0) {
        const sphereGeom = new THREE.SphereGeometry(0.6, 16, 16); // Slightly larger
        const sphereMat = new THREE.MeshBasicMaterial({ color: 0x3b82f6 });
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                const sp = new THREE.Mesh(sphereGeom, sphereMat.clone());
                sp.userData = { i, j };
                pointMeshes.push(sp);
            }
        }
    }

    pointMeshes.forEach(sp => {
        const { i, j } = sp.userData;
        const cp = patch.controlPoints[i][j];
        // Only update pos if it's NOT currently being manipulated
        if (!transformControls.object || transformControls.object !== sp) {
            sp.position.set(cp.x, cp.y, cp.z);
        }
        gridGroup.add(sp);
    });

    // Control Lattice Lines
    for (let i = 0; i < 4; i++) {
        const ptsU = [], ptsV = [];
        for (let j = 0; j < 4; j++) {
            const cpU = patch.controlPoints[i][j];
            const cpV = patch.controlPoints[j][i];
            ptsU.push(cpU.x, cpU.y, cpU.z); ptsV.push(cpV.x, cpV.y, cpV.z);
        }
        gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(ptsU, 3)), lineMat));
        gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(ptsV, 3)), lineMat));
    }
}

// --- Interaction ---
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

function updateMouse(event) {
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
}

function onPointerDown(event) {
    updateMouse(event);
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(pointMeshes);
    
    if (intersects.length > 0) {
        const selected = intersects[0].object;
        console.log("Selected point:", selected.userData);
        transformControls.attach(selected);
        pointMeshes.forEach(p => p.material.color.set(0x3b82f6)); 
        selected.material.color.set(0xffffff); // Active white
    } else if (!transformControls.dragging) {
        transformControls.detach();
        pointMeshes.forEach(p => p.material.color.set(0x3b82f6));
    }
}

function onMouseMove(event) {
    updateMouse(event);
    if (transformControls.dragging) {
        renderer.domElement.style.cursor = 'grabbing';
        return;
    }
    
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(pointMeshes);
    
    // Hover effect
    let isHovering = false;
    pointMeshes.forEach(p => {
        if (transformControls.object !== p) {
            p.material.color.set(0x3b82f6);
            p.scale.set(1, 1, 1);
        }
    });

    if (intersects.length > 0) {
        const hovered = intersects[0].object;
        if (transformControls.object !== hovered) {
            hovered.material.color.set(0x60a5fa);
            hovered.scale.set(1.4, 1.4, 1.4); // More obvious hover
            renderer.domElement.style.cursor = 'pointer';
            isHovering = true;
        }
    } else {
        renderer.domElement.style.cursor = 'default';
    }
}

transformControls.addEventListener('dragging-changed', (e) => {
    orbitControls.enabled = !e.value;
});

transformControls.addEventListener('objectChange', () => {
    const obj = transformControls.object;
    if (obj) {
        const { i, j } = obj.userData;
        patch.controlPoints[i][j] = { x: obj.position.x, y: obj.position.y, z: obj.position.z };
        createSurface();
        updateGrid();
    }
});

document.getElementById('reset-btn').onclick = () => {
    initDefaultPoints();
    pointMeshes = []; 
    transformControls.detach();
    createSurface();
    updateGrid();
};

// --- Standard Events ---
document.getElementById('resolution-slider').oninput = createSurface;
document.getElementById('show-grid').onchange = updateGrid;
document.getElementById('show-wireframe').onchange = createSurface;

// Use Pointer events for better cross-platform support
renderer.domElement.addEventListener('pointerdown', onPointerDown);
renderer.domElement.addEventListener('pointermove', onMouseMove);

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

// Start
createSurface();
updateGrid();
animate();
