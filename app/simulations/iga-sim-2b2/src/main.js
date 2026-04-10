/**
 * Phase 2B.2 | Linear IGA Benchmark
 * Infinite Plate with a Circular Hole
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

// Benchmark Constants
const R = 1.0;
const L = 4.0;
solver.E = 100000;
solver.nu = 0.3;

// The preset now generates the correct Top-Right quadrant directly
let basePatch = NURBSPresets.generatePlateWithHole(R, L);
let patch = { ...basePatch };

let analysisData = {
    u: null,
    load: 100
};

let targetState = { 
    load: 100,
    h: 0,
    p: 2,
    k: false,
    E: 100000,
    nu: 0.30,
    useHybridBasis: true
};
let activeState = { load: -1, h: -1, p: -1, k: false, E: -1, nu: -1 };

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
    
    // Show warning if complexity is high
    const warning = document.getElementById('refinement-warning');
    if (warning) {
        warning.classList.toggle('hidden', targetState.h < 3 && targetState.p < 4);
    }
}

function applyRefinements() {
    // Reset to Base
    patch = JSON.parse(JSON.stringify(basePatch));
    
    // Apply Refinement Stack (p elevation then h subdivision)
    RefineUtils.apply(engine, patch, { 
        h: targetState.h, 
        p: targetState.p 
    });
    
    updateStatusLabels();
    updateSurface();
    updateBoundaryVisuals();
    updateControlPoints();
    updateForceArrows();
}

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

function precomputeReferenceCache(u_ref) {
    const t0 = performance.now();
    const res = 50;
    const nPts = (res + 1) * (res + 1);
    const stress_ref = new Float32Array(nPts);
    const dx_ref = new Float32Array(nPts);
    const dy_ref = new Float32Array(nPts);
    const px_ref = new Float32Array(nPts);
    const py_ref = new Float32Array(nPts);
    const pz_ref = new Float32Array(nPts);

    const stressCap = 4.0 * Math.abs(solverCache.referenceData.load_ref);

    for (let i = 0; i <= res; i++) {
        for (let j = 0; j <= res; j++) {
            const u = i / res;
            const v = j / res;
            const idx = i * (res + 1) + j;
            
            const state = engine.getSurfaceState(patch, u, v);
            px_ref[idx] = state.position.x;
            py_ref[idx] = state.position.y;
            pz_ref[idx] = state.position.z;

            // Displacement
            const interp = interpolateDisplacement(u, v, state.denominator, u_ref);
            if (Number.isFinite(interp.x) && Number.isFinite(interp.y)) {
                dx_ref[idx] = interp.x;
                dy_ref[idx] = interp.y;
            } else {
                dx_ref[idx] = 0;
                dy_ref[idx] = 0;
            }

            const s = solver.getNumericalStress(patch, u_ref, u, v, targetState.E, targetState.nu);
            const val = Number.isFinite(s.vonMises) ? Math.min(s.vonMises, stressCap) : 0;
            stress_ref[idx] = val;
        }
    }
    solverCache.referenceData.stress_ref = stress_ref;
    solverCache.referenceData.dx_ref = dx_ref;
    solverCache.referenceData.dy_ref = dy_ref;
    solverCache.referenceData.px_ref = px_ref;
    solverCache.referenceData.py_ref = py_ref;
    solverCache.referenceData.pz_ref = pz_ref;
    console.log(`Vis Cache Rebuilt in ${(performance.now() - t0).toFixed(1)}ms`);
}

function updateSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 50;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    const mode = viewMode;
    let maxVal = 1;

    // SCALAR SCALING FACTOR
    const currentLoad = targetState.load;
    const ref = solverCache.referenceData;
    const scale = (ref && solverCache.useScalarScaling) ? (currentLoad / ref.load_ref) : 1.0;
    const u_active = (ref && solverCache.useScalarScaling) ? ref.u_ref : (analysisData.u || null);

    if (u_active) {
        if (mode === 'displacement') maxVal = calculateMaxDisp(u_active) * Math.abs(scale);
        else maxVal = calculateMaxStress(u_active) * Math.abs(scale);
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

            if (u_active && ref && ref.dx_ref) {
                // EXTREMELY FAST CACHE LOOKUP
                const idx = i * (res + 1) + j;
                const dx = ref.dx_ref[idx] * scale;
                const dy = ref.dy_ref[idx] * scale;
                
                const px = ref.px_ref[idx] + dx;
                const py = ref.py_ref[idx] + dy;
                const pz = ref.pz_ref[idx];
                positions.push(px, py, pz);

                if (mode === 'displacement') {
                    val = Math.abs(dx); 
                } else if (ref.stress_ref) {
                    val = ref.stress_ref[idx] * Math.abs(scale);
                }
                const t = Math.min(Math.max(val / (maxVal || 1e-6), 0), 1.0);
                const color = getHeatmapColor(t);
                colors.push(color.r, color.g, color.b);
            } else {
                // FALLBACK
                const state = engine.getSurfaceState(patch, u, v);
                let p = { ...state.position };

                if (u_active) {
                    const interp = interpolateDisplacement(u, v, state.denominator, u_active);
                    const dx = interp.x * scale;
                    const dy = interp.y * scale;

                    if (Number.isFinite(dx) && Number.isFinite(dy)) {
                        p.x += dx;
                        p.y += dy;
                    }
                    
                    if (mode === 'displacement') {
                        val = Math.abs(dx); 
                    } else if (ref && ref.stress_ref) {
                        val = ref.stress_ref[i * (res + 1) + j] * Math.abs(scale);
                    } else {
                        const sRef = solver.getNumericalStress(patch, u_active, u, v, targetState.E, targetState.nu);
                        const stressCap = 4.0 * Math.abs(ref ? ref.load_ref : currentLoad);
                        val = Number.isFinite(sRef.vonMises) ? (Math.min(sRef.vonMises, stressCap) * Math.abs(scale)) : 0;
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
    
    // Bottom edge (Angular i=0 is x-axis) -> Fix Y
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[0][j];
        const marker = createBCMarker('y');
        marker.position.set(cp.x, cp.y, cp.z);
        scene.add(marker);
        boundaryVisuals.push(marker);
    }
    // Left edge (Angular i=nU-1 is y-axis) -> Fix X
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[nU-1][j];
        const marker = createBCMarker('x');
        marker.position.set(cp.x, cp.y, cp.z);
        scene.add(marker);
        boundaryVisuals.push(marker);
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

    // Dynamically detect the actual physical extent of the boundary
    const posStart = engine.evaluateSurface(patch, 0.0, 1.0);
    const posEnd = engine.evaluateSurface(patch, 0.5, 1.0); // u=0.5 is the top corner
    const yMin = Math.min(posStart.y, posEnd.y);
    const yMax = Math.max(posStart.y, posEnd.y);

    // Use spatial sampling for uniform visualization (Linear steps in actual physical space)
    const nArrows = 8;
    for (let i = 0; i <= nArrows; i++) {
        const yTarget = yMin + (i / nArrows) * (yMax - yMin);
        
        // Find u corresponding to yTarget along v=1.0 via binary search or simple dense map
        let bestU = 0;
        let minDist = Infinity;
        const resolution = 200;
        for (let j = 0; j <= resolution; j++) {
            const uTry = j / resolution;
            const pos = engine.evaluateSurface(patch, uTry, 1.0);
            const dist = Math.abs(pos.y - yTarget);
            if (Math.abs(pos.x - L) < 0.2 && dist < minDist) {
                minDist = dist;
                bestU = uTry;
            }
        }

        const state = engine.getSurfaceState(patch, bestU, 1.0);
        const cp = state.position;
        const dir = new THREE.Vector3(1, 0, 0);
        const arrow = new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(cp.x, cp.y, cp.z), targetState.load/5000 + 0.3, 0xef4444);
        scene.add(arrow);
        forceVisuals.push(arrow);
    }
}

async function solverLoop() {
    const currentMeshKey = getMeshKey();
    const meshChanged = targetState.h !== activeState.h || 
                         targetState.p !== activeState.p ||
                         targetState.k !== activeState.k ||
                         targetState.E !== activeState.E ||
                         targetState.nu !== activeState.nu;
    
    const loadChanged = targetState.load !== activeState.load;
    
    // Optimization: If Scalar Scaling is on and only load changed, just re-render
    if (solverCache.useScalarScaling && !meshChanged && loadChanged) {
        updateSurface();
        activeState.load = targetState.load;
    }

    const stateChanged = meshChanged || loadChanged;

    if (!isSolving && stateChanged) {
        isSolving = true;

        if (meshChanged) {
            applyRefinements();
            solverCache.meshKey = currentMeshKey;
        }

        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        const bcs = [];
        if (patch.controlPoints[0]) for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
        if (patch.controlPoints[nU-1]) for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });

        // Solve for Reference Load (or current if Scaling is off)
        const solverLoad = solverCache.useScalarScaling ? 100 : targetState.load;
        const integratedForces = solver.calculateNodalTraction(patch, solverLoad, 'right');
        const loads = [];
        for (let i = 0; i < integratedForces.length; i++) {
            if (integratedForces[i] !== 0) {
                const nodeIdx = Math.floor(i / 2);
                const a = Math.floor(nodeIdx / nV), b = nodeIdx % nV;
                if (i % 2 === 0) loads.push({ i: a, j: b, fx: integratedForces[i], fy: 0 });
                else            loads.push({ i: a, j: b, fx: 0, fy: integratedForces[i] });
            }
        }

        try {
            // Check LU Cache
            let u_sol;
            const useLU = document.getElementById('toggle-lu')?.checked;
            if (useLU) {
                if (!solverCache.LU || meshChanged) {
                    if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = "Factorizing Stiffness...";
                    const K_full = solver.assembleStiffness(patch);
                    solver.applyPenaltyConstraints(K_full, patch);
                    const nDofs = nU * nV * 2;
                    const freeIndices = [], fixedValues = new Map();
                    bcs.forEach(bc => {
                        const baseIdx = (bc.i * nV + bc.j) * 2;
                        if (bc.axis === 'x' || bc.axis === 'both') fixedValues.set(baseIdx, bc.value);
                        if (bc.axis === 'y' || bc.axis === 'both') fixedValues.set(baseIdx + 1, bc.value);
                    });
                    for (let i = 0; i < nDofs; i++) if (!fixedValues.has(i)) freeIndices.push(i);
                    const K_red = Array.from({ length: freeIndices.length }, () => new Float64Array(freeIndices.length));
                    for (let i = 0; i < freeIndices.length; i++) for (let j = 0; j < freeIndices.length; j++) K_red[i][j] = K_full[freeIndices[i]][freeIndices[j]];
                    solverCache.LU = K_red;
                    solverCache.P = solver.decomposeLU(solverCache.LU);
                    solverCache.freeIndices = freeIndices;
                    solverCache.fixedValues = fixedValues;
                }
                const F_full = new Float64Array(nU * nV * 2).fill(0);
                for (let i = 0; i < integratedForces.length; i++) F_full[i] = integratedForces[i];
                const F_red = new Float64Array(solverCache.freeIndices.length);
                for (let i = 0; i < solverCache.freeIndices.length; i++) F_red[i] = F_full[solverCache.freeIndices[i]];
                const sol_red = solver.solveLU(solverCache.LU, solverCache.P, F_red);
                u_sol = new Float64Array(nU * nV * 2);
                let freeIdx = 0;
                for (let i = 0; i < u_sol.length; i++) u_sol[i] = solverCache.fixedValues.has(i) ? solverCache.fixedValues.get(i) : sol_red[freeIdx++];
            } else {
                solver.E = targetState.E; solver.nu = targetState.nu;
                u_sol = solver.solve(patch, bcs, loads);
            }

            if (solverCache.useScalarScaling) {
                solverCache.referenceData = { u_ref: u_sol, load_ref: solverLoad };
                precomputeReferenceCache(u_sol);
                analysisData.u = null; 
            } else {
                analysisData.u = u_sol;
            }

            const t1 = performance.now();
            const u_active = solverCache.useScalarScaling ? solverCache.referenceData.u_ref : u_sol;
            const activeLoad = solverCache.useScalarScaling ? solverCache.referenceData.load_ref : targetState.load;
            const scale = targetState.load / activeLoad;

            // EXCLUSIVE CALCULATION & UI UPDATE
            if (viewMode === 'displacement') {
                const maxD = calculateMaxDisp(u_active) * scale;
                const l2ErrorDisp = solver.calculateRelativeL2DisplacementError(patch, u_active, targetState.E, targetState.nu, activeLoad, 1.0);
                if (document.getElementById('l2-error')) document.getElementById('l2-error').textContent = `${(l2ErrorDisp * 100).toFixed(3)}%`;
                if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `[LU+Scaling] L2 Err = ${(l2ErrorDisp*100).toFixed(4)}% | ${(t1-t0).toFixed(1)}ms`;
            } else {
                const l2Error = solver.calculateRelativeL2Error(patch, u_active, targetState.E, targetState.nu, activeLoad, 1.0);
                if (document.getElementById('l2-error')) document.getElementById('l2-error').textContent = `${(l2Error * 100).toFixed(3)}%`;
                if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `[LU+Scaling] Stress Err = ${(l2Error*100).toFixed(4)}% | ${(t1-t0).toFixed(1)}ms`;
            }

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
    // Load Slider
    document.getElementById('load-slider').oninput = (e) => {
        const val = parseFloat(e.target.value);
        if (document.getElementById('load-val')) document.getElementById('load-val').textContent = `${val} N/mm`;
        targetState.load = val;
    };

    // Material Sliders
    if (document.getElementById('e-slider')) {
        document.getElementById('e-slider').oninput = (e) => {
            const val = parseFloat(e.target.value);
            if (document.getElementById('e-val')) document.getElementById('e-val').textContent = val >= 10000 ? `${(val / 1000).toFixed(0)}k` : val;
            targetState.E = val;
        };
    }
    if (document.getElementById('nu-slider')) {
        document.getElementById('nu-slider').oninput = (e) => {
            const val = parseFloat(e.target.value);
            if (document.getElementById('nu-val')) document.getElementById('nu-val').textContent = val.toFixed(2);
            targetState.nu = val;
        };
    }

    // Refinement Sliders
    if (document.getElementById('h-slider')) {
        document.getElementById('h-slider').oninput = (e) => {
            targetState.h = parseInt(e.target.value);
        };
    }
    if (document.getElementById('p-slider')) {
        document.getElementById('p-slider').oninput = (e) => {
            targetState.p = parseInt(e.target.value);
        };
    }

    // Reset Button
    if (document.getElementById('reset-ref')) {
        document.getElementById('reset-ref').onclick = () => {
            targetState.h = 0;
            targetState.p = 2;
            if (document.getElementById('h-slider')) document.getElementById('h-slider').value = 0;
            if (document.getElementById('p-slider')) document.getElementById('p-slider').value = 2;
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

function calculateMaxStress(u_disp) {
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
            
            const s = solver.getNumericalStress(patch, u_disp, u, v, targetState.E, targetState.nu);
            if (isFinite(s.vonMises) && s.vonMises > max) max = s.vonMises;
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

init();
