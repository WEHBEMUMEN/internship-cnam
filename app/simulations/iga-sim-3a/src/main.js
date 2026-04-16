/**
 * Phase 2B.2 | Linear IGA Benchmark
 * Infinite Plate with a Circular Hole
 */

const engine = new NURBS2D();
const solver = new IGANonlinearSolver(engine);

let currentGeometry = 'cantilever';
function createCantileverBeam(L, H) {
    const p = 2, q = 2;
    const U = [0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1];
    const V = [0, 0, 0, 0.5, 1, 1, 1];
    
    const nu = 6; // 6 control points in U -> gives 6-2=4 elements
    const nv = 4; // 4 control points in V -> gives 4-2=2 elements
    
    const controlPoints = [];
    const weights = [];
    for (let i = 0; i < nu; i++) {
        controlPoints[i] = [];
        weights[i] = [];
        for (let j = 0; j < nv; j++) {
            controlPoints[i][j] = {
                x: (i / (nu - 1)) * L,
                y: (j / (nv - 1)) * H,
                z: 0
            };
            weights[i][j] = 1;
        }
    }
    return { p, q, U, V, controlPoints, weights };
}

let basePatch = createCantileverBeam(10, 2); 
let patch = { ...basePatch };

let analysisData = {
    u: null,
    load: 500
};

let targetState = { 
    load: 500,
    steps: 5,
    iters: 10,
    type: 'linear',
    E: 200000, // Steel
    nu: 0.30,
    geometry: 'cantilever'
};
let activeState = { load: -1, steps: -1, iters: -1, type: '', geometry: '' };

// Solver & Scalar Scaling State
let solverCache = {
    meshKey: "",
    LU: null,
    P: null,
    referenceData: null, // { u_ref: Float64Array, load_ref: number, stress_ref: Float32Array }
    useScalarScaling: true 
};

function getMeshKey() {
    return `${targetState.h}-${targetState.p}-${targetState.k}-${targetState.E}-${targetState.nu}`;
}

let isSolving = false;

function updateStatusLabels() {
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    if (document.getElementById('status-degree')) document.getElementById('status-degree').textContent = `p=${patch.p}, q=${patch.q}`;
    if (document.getElementById('status-spans')) document.getElementById('status-spans').textContent = `${patch.U.length - 2 * patch.p - 1} x ${patch.V.length - 2 * patch.q - 1}`;
    if (document.getElementById('status-points')) document.getElementById('status-points').textContent = `${nU} x ${nV}`;
}

function updateGeometry() {
    if (targetState.geometry === 'cantilever') {
        basePatch = createCantileverBeam(10, 2);
    } else {
        basePatch = NURBSPresets.generatePlateWithHole(1.0, 4.0);
    }
    patch = JSON.parse(JSON.stringify(basePatch));
    updateStatusLabels();
    updateSurface();
    updateBoundaryVisuals();
    updateControlPoints();
    updateForceArrows();
}

// No logic here

function init() {
    setupUI();
    updateSurface();
    updateBoundaryVisuals();
    updateForceArrows();
    animate();
}

// --- Three.js ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);

// Top-Down Orthographic-like Perspective
const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(2.5, 2.5, 12); 
camera.lookAt(2.5, 2.5, 0);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableRotate = false; // Stick to 2D
controls.target.set(2.5, 2.5, 0);

scene.add(new THREE.AmbientLight(0xffffff, 1.0));

// XY orientation helper
const axesHelper = new THREE.AxesHelper(5);
scene.add(axesHelper);

let surfaceMesh, wireframeOverlay;
let boundaryVisuals = [];
let controlPointVisuals = [];
let forceVisuals = [];

let visibilityState = {
    cp: false,
    bc: true,
    force: true
};

let viewMode = 'displacement'; // 'displacement' or 'stress'

window.setViewMode = (m) => {
    viewMode = m;
    
    // Update UI Legend labels
    const modeTitle = document.getElementById('legend-mode-title');
    const modeSub = document.getElementById('legend-mode-sub');
    if (modeTitle) modeTitle.innerText = m === 'displacement' ? 'Displacement' : 'Stress';
    if (modeSub) modeSub.innerText = m === 'displacement' ? 'X-Magnitude' : 'von Mises (MPa)';

    updateSurface();
};

function updateSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 50;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    const mode = viewMode;
    let maxVal = 1;

    const u_active = analysisData.u;

    if (u_active) {
        const isLinear = targetState.type === 'linear';
        const maxD = calculateMaxDisp(u_active);
        const maxS = calculateMaxStress(u_active, isLinear);
        if (mode === 'displacement') maxVal = maxD;
        else maxVal = maxS;

        if (document.getElementById('max-disp')) document.getElementById('max-disp').innerText = `${maxD.toFixed(3)} mm`;
        if (document.getElementById('max-stress')) document.getElementById('max-stress').innerText = `${maxS.toFixed(2)} MPa`;
        
        if (document.getElementById('l2-error')) {
            if (analysisData.l2Diff !== undefined) {
                document.getElementById('l2-error').innerText = `${(analysisData.l2Diff * 100).toFixed(2)}%`;
            } else {
                document.getElementById('l2-error').innerText = "0.00%";
            }
        }
    }

    // Dynamic UI Legend Synchronization
    const legendMax = document.getElementById('legend-val-max');
    const legendMin = document.getElementById('legend-val-min');
    const unit = mode === 'displacement' ? 'mm' : 'MPa';
    if (legendMax) legendMax.innerText = `${maxVal.toFixed(mode === 'displacement' ? 3 : 1)} ${unit}`;
    if (legendMin) legendMin.innerText = `${(0).toFixed(mode === 'displacement' ? 3 : 1)} ${unit}`;

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            let val = 0;

            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };

            if (u_active) {
                const interp = interpolateDisplacement(u, v, state.denominator, u_active);
                const dx = interp.x;
                const dy = interp.y;

                if (Number.isFinite(dx) && Number.isFinite(dy)) {
                    p.x += dx;
                    p.y += dy;
                }
                
                if (mode === 'displacement') {
                    val = Math.abs(dx); 
                } else {
                    const sRef = solver.getNumericalStress(patch, u_active, u, v, targetState.E, targetState.nu, targetState.type === 'linear');
                    const stressCap = 4.0 * targetState.load;
                    val = Number.isFinite(sRef.vonMises) ? Math.min(sRef.vonMises, stressCap) : 0;
                }
                
                const t = Math.min(Math.max(val / (maxVal || 1e-6), 0), 1.0);
                const color = getHeatmapColor(t);
                colors.push(color.r, color.g, color.b);
            } else {
                colors.push(0.1, 0.2, 0.4);
            }
            positions.push(p.x, p.y, p.z);
        }
    }

    const indices = [];
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            const a = i * (res+1) + j;
            const b = (i+1) * (res+1) + j;
            const c = (i+1) * (res+1) + (j+1);
            const d = i * (res+1) + (j+1);
            indices.push(a, b, d, b, c, d);
        }
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

    const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);

    wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.15 })
    );
    scene.add(wireframeOverlay);
}

function interpolateDisplacement(u, v, denom, u_disp) {
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    let dx = 0, dy = 0;
    for (let i = 0; i < nU; i++) {
        const Ni = engine.basis1D(i, patch.p, patch.U, u);
        if (Ni === 0) continue;
        for (let j = 0; j < nV; j++) {
            const Mj = engine.basis1D(j, patch.q, patch.V, v);
            const R = (Ni * Mj * patch.weights[i][j]) / denom;
            const idx = (i * nV + j) * 2;
            dx += R * u_disp[idx];
            dy += R * u_disp[idx + 1];
        }
    }
    return { x: dx, y: dy };
}

function calculateMaxDisp(u) {
    if (!u) return 1;
    let max = 0;
    for (let i = 0; i < u.length; i += 2) {
        const m = Math.abs(u[i]); // X-displacement magnitude as it was before
        if (m > max) max = m;
    }
    return max;
}

function calculateMaxStress(u_disp) {
    if (!u_disp) return 1;
    let max = 0;
    const res = 15; // Higher sampling for peak capture
    for (let i = 0; i <= res; i++) {
        for (let j = 0; j <= res; j++) {
            // Sample angularly (u) 0 to 1, but radially (v) 0 to 0.95
            // The degenerate corner at v=1 is a numerical artifact, the physics is at v=0 (hole).
            const u = i / res;
            const v = (j / res) * 0.95; 
            const s = solver.getNumericalStress(patch, u_disp, u, v, targetState.E, targetState.nu);
            if (s.vonMises > max) max = s.vonMises;
        }
    }
    return max || 1;
}

function getHeatmapColor(t) {
    // Standard Jet/Rainbow scale matching the UI legend: Blue -> Cyan -> Green -> Yellow -> Orange -> Red
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

function createBCMarker(type = 'x') {
    const group = new THREE.Group();
    const triGeo = new THREE.ConeGeometry(0.12, 0.2, 3);
    const triMat = new THREE.MeshBasicMaterial({ color: 0x334155 });
    const triangle = new THREE.Mesh(triGeo, triMat);
    const rollGeo = new THREE.CircleGeometry(0.04, 16);
    const rollMat = new THREE.MeshBasicMaterial({ color: 0x64748b });
    const roller = new THREE.Mesh(rollGeo, rollMat);
    
    if (type === 'both') {
        const wallGeo = new THREE.PlaneGeometry(0.2, 0.5);
        const wallMat = new THREE.MeshBasicMaterial({ color: 0x64748b });
        const wall = new THREE.Mesh(wallGeo, wallMat);
        wall.position.x = -0.1;
        group.add(wall);
        return group;
    }
    
    if (type === 'y') {
        triangle.position.y = -0.15;
        roller.position.y = -0.28;
    } else {
        triangle.rotation.z = -Math.PI / 2; // Inverted to point right -> left
        triangle.position.x = -0.15;
        roller.position.x = -0.28;
    }
    
    group.add(triangle); group.add(roller);
    return group;
}

function updateBoundaryVisuals() {
    boundaryVisuals.forEach(obj => scene.remove(obj));
    boundaryVisuals = [];
    if (!visibilityState.bc) return;

    const nU = patch.controlPoints.length; 
    const nV = patch.controlPoints[0].length; 
    
    if (targetState.geometry === 'cantilever') {
        // Clamp Left edge (i=0)
        for(let j=0; j<nV; j++) {
            const cp = patch.controlPoints[0][j];
            const marker = createBCMarker('both');
            marker.position.set(cp.x, cp.y, cp.z);
            scene.add(marker);
            boundaryVisuals.push(marker);
        }
    } else {
        // Plate with Hole symmetry
        // Bottom edge (i=0) -> Fix Y
        for(let j=0; j<nV; j++) {
            const cp = patch.controlPoints[0][j];
            const marker = createBCMarker('y');
            marker.position.set(cp.x, cp.y, cp.z);
            scene.add(marker);
            boundaryVisuals.push(marker);
        }
        // Left edge (i=nU-1 is the vertical symmetry line) -> Fix X
        for(let j=0; j<nV; j++) {
            const cp = patch.controlPoints[nU-1][j];
            const marker = createBCMarker('x');
            marker.position.set(cp.x, cp.y, cp.z);
            scene.add(marker);
            boundaryVisuals.push(marker);
        }
    }
}

function updateControlPoints() {
    controlPointVisuals.forEach(obj => scene.remove(obj));
    controlPointVisuals = [];
    if (!visibilityState.cp) return;

    const geo = new THREE.SphereGeometry(0.06);
    const mat = new THREE.MeshBasicMaterial({ color: 0x10b981 });
    patch.controlPoints.forEach((row, i) => {
        row.forEach((cp, j) => {
            const sphere = new THREE.Mesh(geo, mat);
            sphere.position.set(cp.x, cp.y, cp.z);
            scene.add(sphere);
            controlPointVisuals.push(sphere);
        });
    });
}

function updateForceArrows() {
    forceVisuals.forEach(obj => scene.remove(obj));
    forceVisuals = [];
    if (!visibilityState.force || targetState.load === 0) return;

    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;

    let visualLoads = [];
    if (targetState.geometry === 'cantilever') {
        const p1 = patch.controlPoints[nU - 1][0];
        const p2 = patch.controlPoints[nU - 1][nV - 1];
        visualLoads.push({ pos: p1, dir: new THREE.Vector3(0, -1, 0) });
        visualLoads.push({ pos: p2, dir: new THREE.Vector3(0, -1, 0) });
    } else {
        const p1 = patch.controlPoints[0][nV - 1];
        const p2 = patch.controlPoints[1][nV - 1];
        visualLoads.push({ pos: p1, dir: new THREE.Vector3(1, 0, 0) });
        visualLoads.push({ pos: p2, dir: new THREE.Vector3(1, 0, 0) });
    }

    visualLoads.forEach(load => {
        const arrow = new THREE.ArrowHelper(load.dir, new THREE.Vector3(load.pos.x, load.pos.y, load.pos.z), 1.0, 0xef4444);
        scene.add(arrow);
        forceVisuals.push(arrow);
    });
}

async function solverLoop() {
    const stateChanged = 
        targetState.type !== activeState.type || 
        targetState.steps !== activeState.steps ||
        targetState.iters !== activeState.iters ||
        targetState.geometry !== activeState.geometry;
    
    const loadChanged = targetState.load !== activeState.load;

    if (!isSolving && (stateChanged || loadChanged)) {
        isSolving = true;

        if (targetState.geometry !== activeState.geometry) {
            updateGeometry();
        }

        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        const bcs = [];
        if (targetState.geometry === 'cantilever') {
            // Fix left edge
            for (let j = 0; j < nV; j++) {
                bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
            }
        } else {
            // Symmetry for plate with hole
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
            for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j: j, axis: 'x', value: 0 });
        }

        const loads = [];
        if (targetState.geometry === 'cantilever') {
            // Load on right edge, tip
            loads.push({ type: 'nodal', i: nU - 1, j: nV - 1, fx: 0, fy: -targetState.load });
            loads.push({ type: 'nodal', i: nU - 1, j: 0, fx: 0, fy: -targetState.load });
        } else {
            // Uniaxial tension on the right edge corresponding to X=L (i=0, i=1)
            loads.push({ type: 'nodal', i: 0, j: nV - 1, fx: targetState.load, fy: 0 });
            loads.push({ type: 'nodal', i: 1, j: nV - 1, fx: targetState.load, fy: 0 });
        }

        try {
            solver.E = targetState.E;
            solver.nu = targetState.nu;

            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Running ${targetState.type} Analysis...`;

            const res = solver.solveNonlinear(patch, bcs, loads, {
                iterations: targetState.iters,
                tolerance: 1e-5,
                steps: targetState.type === 'linear' ? 1 : targetState.steps,
                isLinear: targetState.type === 'linear',
                onProgress: (info) => {
                    if (document.getElementById('conv-stat')) {
                        document.getElementById('conv-stat').textContent = `Step ${info.step}/${targetState.steps} | Iter ${info.iter} | Res: ${info.norm.toExponential(2)}`;
                    }
                }
            });

            let l2Diff = 0;
            if (targetState.type !== 'linear') {
                const linearRes = solver.solveNonlinear(patch, bcs, loads, { iterations: 1, steps: 1, isLinear: true });
                let sumDiff = 0, sumNorm = 0;
                for(let i=0; i<linearRes.u.length; i++) {
                    const d = res.u[i] - linearRes.u[i];
                    sumDiff += d*d;
                    sumNorm += linearRes.u[i]*linearRes.u[i];
                }
                l2Diff = sumNorm > 0 ? Math.sqrt(sumDiff) / Math.sqrt(sumNorm) : 0;
            }

            analysisData.u = res.u;
            analysisData.l2Diff = l2Diff;
            
            const t1 = performance.now();
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Analysis Complete in ${(t1 - t0).toFixed(1)}ms. Max Residual: ${res.residualHistory[res.residualHistory.length - 1]?.norm.toExponential(2) || '0'}`;

            updateSurface();
            updateBoundaryVisuals();
            updateControlPoints();
            updateForceArrows();
            
            activeState = JSON.parse(JSON.stringify(targetState));
        } catch (err) {
            console.error(err);
        }
        isSolving = false;
    }
    requestAnimationFrame(solverLoop);
}

function setupUI() {
    window.setAnalysisTypeHandler = (t) => { targetState.type = t; };

    // Load Slider
    document.getElementById('load-slider').oninput = (e) => {
        const val = parseFloat(e.target.value);
        if (document.getElementById('load-val')) document.getElementById('load-val').textContent = `${val}`;
    };
    
    // Geometry Select
    document.getElementById('geometry-select').onchange = (e) => {
        targetState.geometry = e.target.value;
    };
    
    document.getElementById('steps-slider').oninput = (e) => {
        const val = parseInt(e.target.value);
        if (document.getElementById('steps-val')) document.getElementById('steps-val').textContent = val;
        targetState.steps = val;
    };
    
    document.getElementById('iters-slider').oninput = (e) => {
        const val = parseInt(e.target.value);
        if (document.getElementById('iters-val')) document.getElementById('iters-val').textContent = val;
        targetState.iters = val;
    };

    // Calculate Button
    if (document.getElementById('solve-btn')) {
        document.getElementById('solve-btn').onclick = () => {
            const loadEl = document.getElementById('load-slider');
            if (loadEl) targetState.load = parseFloat(loadEl.value);
        };
    }

    // Reset Button
    if (document.getElementById('reset-ref')) {
        document.getElementById('reset-ref').onclick = () => {
            targetState.load = 0;
            if (document.getElementById('load-slider')) document.getElementById('load-slider').value = 0;
        };
    }

    // K-Refinement Shortcut
    if (document.getElementById('apply-k')) {
        document.getElementById('apply-k').onclick = () => {
            targetState.p = Math.min(targetState.p + 1, 5);
            targetState.h = Math.min(targetState.h + 1, 4);
            if (document.getElementById('p-slider')) document.getElementById('p-slider').value = targetState.p;
            if (document.getElementById('h-slider')) document.getElementById('h-slider').value = targetState.h;
        };
    }
    
    if (document.getElementById('generate-plot-btn')) {
        document.getElementById('generate-plot-btn').onclick = generateFDPlot;
    }
    
    document.getElementById('toggle-cp').onchange = (e) => {
        visibilityState.cp = e.target.checked;
        updateControlPoints();
    };
    document.getElementById('toggle-bc').onchange = (e) => {
        visibilityState.bc = e.target.checked;
        updateBoundaryVisuals();
    };
    document.getElementById('toggle-force').onchange = (e) => {
        visibilityState.force = e.target.checked;
        updateForceArrows();
    };

    updateStatusLabels();
    solverLoop();
}

function calculateMaxDisp(u) {
    let max = 0;
    for (let i = 0; i < u.length; i+=2) {
        const mag = Math.sqrt(u[i]*u[i] + u[i+1]*u[i+1]);
        if (mag > max) max = mag;
    }
    return max;
}

function calculateMaxStress(u_disp, isLinear = false) {
    let max = 0;
    const res = 15;
    const L = 4.0;
    for (let i = 0; i <= res; i++) {
        for (let j = 0; j <= res; j++) {
            const u = i / res;
            const v = j / res;
            // Skip the outer boundary region (degenerate mapping zone)
            const pos = engine.evaluateSurface(patch, u, v);
            if (pos.x > 0.95 * L || pos.y > 0.95 * L) continue;
            
            const s = solver.getNumericalStress(patch, u_disp, u, v, targetState.E, targetState.nu, isLinear);
            if (isFinite(s.vonMises) && s.vonMises > max) max = s.vonMises;
        }
    }
    return max;
}

function calculateMaxStrain(u_disp, isLinear = false) {
    let max = 0;
    const res = 15;
    const L = 4.0;
    for (let i = 0; i <= res; i++) {
        for (let j = 0; j <= res; j++) {
            const u = i / res;
            const v = j / res;
            const pos = engine.evaluateSurface(patch, u, v);
            if (pos.x > 0.95 * L || pos.y > 0.95 * L) continue;
            
            const s = solver.getNumericalStress(patch, u_disp, u, v, targetState.E, targetState.nu, isLinear);
            if (isFinite(s.vonMisesStrain) && s.vonMisesStrain > max) max = s.vonMisesStrain;
        }
    }
    return max;
}

function calculateHMax(patch) {
    const { U, V } = patch;
    const uniqueU = [...new Set(U)];
    const uniqueV = [...new Set(V)];
    let hMax = 0;
    for (let i = 0; i < uniqueU.length - 1; i++) {
        for (let j = 0; j < uniqueV.length - 1; j++) {
            const uMin = uniqueU[i]; const uMax = uniqueU[i+1];
            const vMin = uniqueV[j]; const vMax = uniqueV[j+1];
            if (uMax - uMin < 1e-9 || vMax - vMin < 1e-9) continue;
            
            const p00 = engine.evaluateSurface(patch, uMin, vMin);
            const p11 = engine.evaluateSurface(patch, uMax, vMax);
            const diag = Math.sqrt((p11.x - p00.x)**2 + (p11.y - p00.y)**2);
            if (diag > hMax) hMax = diag;
        }
    }
    return hMax;
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}

window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

let chartInstance = null;
let ssChartInstance = null;

async function generateFDPlot() {
    isSolving = true;
    if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Sweeping Load... 0%`;

    const nPoints = 20;
    const maxLoad = targetState.geometry === 'cantilever' ? 50000 : 100000;
    
    const linearData = [];
    const nonlinearData = [];
    
    const linearSSData = [];
    const nonlinearSSData = [];
    
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    const bcs = [];
    if (targetState.geometry === 'cantilever') {
        for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
    } else {
        if (patch.controlPoints[0]) for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
        if (patch.controlPoints[nU-1]) for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });
    }
    
    const modal = document.getElementById('chart-modal');
    modal.classList.remove('hidden');

    for (let i = 0; i <= nPoints; i++) {
        const load = (i / nPoints) * maxLoad;
        const loads = [];
        if (targetState.geometry === 'cantilever') {
            loads.push({ type: 'nodal', i: nU - 1, j: nV - 1, fx: 0, fy: -load });
            loads.push({ type: 'nodal', i: nU - 1, j: 0, fx: 0, fy: -load });
        } else {
            loads.push({ type: 'nodal', i: 0, j: nV - 1, fx: load, fy: 0 });
            loads.push({ type: 'nodal', i: 1, j: nV - 1, fx: load, fy: 0 });
        }
        
        solver.E = targetState.E;
        solver.nu = targetState.nu;

        // Linear Baseline
        const linRes = solver.solveNonlinear(patch, bcs, loads, { iterations: 1, steps: 1, isLinear: true });
        const linDisp = calculateMaxDisp(linRes.u);
        const linStress = calculateMaxStress(linRes.u, true);
        const linStrain = calculateMaxStrain(linRes.u, true);
        linearData.push({ x: linDisp, y: load });
        linearSSData.push({ x: linStrain, y: linStress });

        // Nonlinear Point
        const nlRes = solver.solveNonlinear(patch, bcs, loads, { iterations: 20, steps: targetState.geometry === 'cantilever' ? 5 : 1, isLinear: false });
        const nlDisp = calculateMaxDisp(nlRes.u);
        const nlStress = calculateMaxStress(nlRes.u, false);
        const nlStrain = calculateMaxStrain(nlRes.u, false);
        nonlinearData.push({ x: nlDisp, y: load });
        nonlinearSSData.push({ x: nlStrain, y: nlStress });
        
        if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Sweeping Load... ${((i/nPoints)*100).toFixed(0)}%`;
        
        // Live update the chart as it solves
        renderChart(linearData, nonlinearData, linearSSData, nonlinearSSData);
        await new Promise(r => setTimeout(r, 10)); // allow UI update
    }
    
    isSolving = false;
    activeState.load = -1; 
}

function renderChart(linearData, nonlinearData, linearSSData, nonlinearSSData) {
    const ctxFD = document.getElementById('fdChart').getContext('2d');
    const ctxSS = document.getElementById('ssChart').getContext('2d');
    
    if (chartInstance) chartInstance.destroy();
    if (ssChartInstance) ssChartInstance.destroy();
    
    chartInstance = new Chart(ctxFD, {
        type: 'line',
        data: {
            datasets: [{
                label: 'Linear Hookean Theory',
                data: linearData,
                borderColor: '#94a3b8',
                borderDash: [5, 5],
                backgroundColor: 'transparent',
                tension: 0.1
            }, {
                label: 'Geometric Nonlinearity',
                data: nonlinearData,
                borderColor: '#10b981',
                backgroundColor: 'rgba(16, 185, 129, 0.1)',
                fill: true,
                tension: 0.4
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            animation: false, // disable animation for live update
            scales: {
                x: {
                    type: 'linear',
                    title: { display: true, text: 'Max Displacement (mm)', color: '#94a3b8' },
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#64748b' }
                },
                y: {
                    type: 'linear', 
                    title: { display: true, text: 'Applied Force', color: '#94a3b8' },
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#64748b' }
                }
            },
            plugins: { 
                legend: { labels: { color: '#cbd5e1' } },
                zoom: {
                    zoom: {
                        wheel: { enabled: true },
                        drag: { enabled: true, backgroundColor: 'rgba(16, 185, 129, 0.2)' },
                        pinch: { enabled: true },
                        mode: 'xy'
                    },
                    pan: {
                        enabled: true,
                        mode: 'xy',
                        modifierKey: 'ctrl'
                    }
                }
            }
        }
    });

    ssChartInstance = new Chart(ctxSS, {
        type: 'line',
        data: {
            datasets: [{
                label: 'Linear Stress Theory',
                data: linearSSData,
                borderColor: '#94a3b8',
                borderDash: [5, 5],
                backgroundColor: 'transparent',
                tension: 0.1
            }, {
                label: 'Nonlinear Von Mises Stress',
                data: nonlinearSSData,
                borderColor: '#f43f5e',
                backgroundColor: 'rgba(244, 63, 94, 0.1)',
                fill: true,
                tension: 0.4
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            animation: false,
            scales: {
                x: {
                    type: 'linear',
                    title: { display: true, text: 'Max Strain (mm/mm)', color: '#94a3b8' },
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#64748b' }
                },
                y: {
                    type: 'linear', 
                    title: { display: true, text: 'Max Von Mises Stress (MPa)', color: '#94a3b8' },
                    grid: { color: 'rgba(255,255,255,0.05)' },
                    ticks: { color: '#64748b' }
                }
            },
            plugins: { 
                legend: { labels: { color: '#cbd5e1' } },
                zoom: {
                    zoom: {
                        wheel: { enabled: true },
                        drag: { enabled: true, backgroundColor: 'rgba(244, 63, 94, 0.2)' },
                        pinch: { enabled: true },
                        mode: 'xy'
                    },
                    pan: {
                        enabled: true,
                        mode: 'xy',
                        modifierKey: 'ctrl'
                    }
                }
            }
        }
    });

    document.getElementById('close-modal-btn').onclick = () => {
        document.getElementById('chart-modal').classList.add('hidden');
    };
    
    document.getElementById('reset-zoom-btn').onclick = () => {
        if (chartInstance) chartInstance.resetZoom();
        if (ssChartInstance) ssChartInstance.resetZoom();
    };
}

init();
