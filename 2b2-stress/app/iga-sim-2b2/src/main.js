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
    useHybridBasis: true,
    useHybridSolver: false
};
let activeState = { load: -1, h: -1, p: -1, k: false, E: -1, nu: -1, useHybridSolver: -1 };
let localNodeOnline = false;

async function checkLocalNode() {
    try {
        const res = await fetch('http://localhost:8000/', { mode: 'cors' });
        const data = await res.json();
        const wasOnline = localNodeOnline;
        localNodeOnline = data.status === 'online';
        
        const statusEl = document.getElementById('node-status');
        const nodeToggle = document.getElementById('toggle-hybrid');
        
        if (statusEl) {
            statusEl.innerText = localNodeOnline ? 'Online' : 'Offline';
            statusEl.className = localNodeOnline ? 'text-emerald-400 font-bold' : 'text-slate-400';
        }
        if (nodeToggle) {
            nodeToggle.disabled = !localNodeOnline;
            if (!localNodeOnline) {
                nodeToggle.checked = false;
                targetState.useHybridSolver = false;
            }
        }
    } catch(e) {
        localNodeOnline = false;
        const statusEl = document.getElementById('node-status');
        if (statusEl) {
            statusEl.innerText = 'Offline';
            statusEl.className = 'text-slate-400';
        }
    }
}
setInterval(checkLocalNode, 2000);
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

function init() {
    setupUI();
    updateSurface();
    updateBoundaryVisuals();
    updateForceArrows();
    animate();
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

    // EXCLUSIVE MAX CALCULATION
    if (analysisData.u) {
        if (mode === 'displacement') maxVal = calculateMaxDisp(analysisData.u);
        else maxVal = calculateMaxStress(analysisData.u);
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
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };

            if (analysisData.u) {
                const interp = interpolateDisplacement(u, v, state.denominator);
                
                // NAN SAFETY FOR GEOMETRY
                if (Number.isFinite(interp.x) && Number.isFinite(interp.y)) {
                    p.x += interp.x;
                    p.y += interp.y;
                }
                
                let val = 0;
                if (mode === 'displacement') {
                    val = Math.abs(interp.x); 
                } else {
                    let s;
                    if (analysisData.remote) {
                        // Use pre-computed remote results if available
                        // Finding nearest point in high-res grid is slow, 
                        // so for "testing" we just interpolate if we had a grid mapper.
                        // For now we'll call local JS for visual continuity unless we strictly want remote values.
                        s = solver.getNumericalStress(patch, analysisData.u, u, v, targetState.E, targetState.nu);
                    } else {
                        s = solver.getNumericalStress(patch, analysisData.u, u, v, targetState.E, targetState.nu);
                    }
                    // Cap stress at a physical maximum to prevent singularity blow-ups
                    // from washing out the heatmap (SCF = 3, plus margin → 4× Tx)
                    const stressCap = 4.0 * Math.abs(targetState.load);
                    val = Number.isFinite(s.vonMises) ? Math.min(s.vonMises, stressCap) : 0;
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

function interpolateDisplacement(u, v, denom) {
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
            dx += R * analysisData.u[idx];
            dy += R * analysisData.u[idx + 1];
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
        const arrow = new THREE.ArrowHelper(dir, new THREE.Vector3(cp.x, cp.y, cp.z), targetState.load/250, 0xef4444);
        scene.add(arrow);
        forceVisuals.push(arrow);
    }
}

async function solverLoop() {
    const stateChanged = targetState.load !== activeState.load || 
                         targetState.h !== activeState.h || 
                         targetState.p !== activeState.p ||
                         targetState.k !== activeState.k ||
                         targetState.E !== activeState.E ||
                         targetState.nu !== activeState.nu;

    if (!isSolving && stateChanged) {
        isSolving = true;

        // Sync Material Properties to Solver
        solver.E = targetState.E;
        solver.nu = targetState.nu;
        
        // If topology changed, re-apply refinements
        if (targetState.h !== activeState.h || targetState.p !== activeState.p || targetState.k !== activeState.k) {
            applyRefinements();
        }

        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        const bcs = [];
        // Edge u=0 is bottom (y=0) -> Fix Y (Symmetry)
        if (patch.controlPoints[0]) {
            for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
        }
        
        // NEW: Physically accurate 1D Gaussian Integration for boundary traction
        const integratedForces = solver.calculateNodalTraction(patch, targetState.load, 'right');
        const loads = [];
        for (let i = 0; i < integratedForces.length; i++) {
            if (integratedForces[i] !== 0) {
                // Map flat index back to {i, j, axis} if needed, or update solver to take raw vector
                const dofId = i % 2;
                const nodeIdx = Math.floor(i / 2);
                const a = Math.floor(nodeIdx / nV);
                const b = nodeIdx % nV;
                
                if (dofId === 0) loads.push({ i: a, j: b, fx: integratedForces[i], fy: 0 });
                else            loads.push({ i: a, j: b, fx: 0, fy: integratedForces[i] });
            }
        }

        // Edge u=1 is left (x=0) -> Fix X (Symmetry)
        if (patch.controlPoints[nU-1]) {
            for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });
        }

        try {
            if (targetState.useHybridSolver && localNodeOnline) {
                // REMOTE PYTHON SOLVE (Nutils)
                const response = await fetch('http://localhost:8000/solve', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        radius: 1.0,
                        length: 4.0,
                        E: targetState.E,
                        nu: targetState.nu,
                        traction: targetState.load,
                        nrefine: targetState.h
                    })
                });
                const remoteData = await response.json();
                // Store remote stress data for visualization
                analysisData.remote = remoteData;
                analysisData.u = new Float64Array(patch.controlPoints.length * patch.controlPoints[0].length * 2); // Dummy for compatibility
            } else {
                // LOCAL JS SOLVE
                analysisData.u = solver.solve(patch, bcs, loads);
                analysisData.remote = null;
            }
            const t1 = performance.now();
            
            // EXCLUSIVE CALCULATION & UI UPDATE
            if (viewMode === 'displacement') {
                const maxD = calculateMaxDisp(analysisData.u);
                if (document.getElementById('max-disp')) document.getElementById('max-disp').textContent = `${maxD.toFixed(4)} mm`;
                if (document.getElementById('max-stress')) document.getElementById('max-stress').textContent = `---`;
                if (document.getElementById('l2-error')) document.getElementById('l2-error').textContent = `---`;
                if (document.getElementById('legend-val-max')) document.getElementById('legend-val-max').innerText = `${maxD.toFixed(3)} mm`;
                if (document.getElementById('legend-val-min')) document.getElementById('legend-val-min').innerText = `0.000 mm`;
            } else {
                // Stress analysis + Benchmark error
                const maxS = calculateMaxStress(analysisData.u);
                const R_hole = 1.0; 
                const l2Error = solver.calculateRelativeL2Error(patch, analysisData.u, targetState.E, targetState.nu, targetState.load, R_hole);

                if (document.getElementById('max-disp')) document.getElementById('max-disp').textContent = `---`;
                if (document.getElementById('max-stress')) document.getElementById('max-stress').textContent = `${maxS.toFixed(2)} MPa`;
                if (document.getElementById('l2-error')) document.getElementById('l2-error').textContent = `${(l2Error * 100).toFixed(3)}%`;
                if (document.getElementById('legend-val-max')) document.getElementById('legend-val-max').innerText = `${maxS.toFixed(1)} MPa`;
                if (document.getElementById('legend-val-min')) document.getElementById('legend-val-min').innerText = `0.0 MPa`;
                
                if (document.getElementById('conv-stat')) {
                    document.getElementById('conv-stat').textContent = `Stress Benchmarked: L2 Err = ${(l2Error*100).toFixed(4)}% | ${(t1 - t0).toFixed(1)}ms`;
                }
            }

            updateSurface();
            updateBoundaryVisuals();
            updateControlPoints();
            updateForceArrows();
            
            if (viewMode === 'displacement' && document.getElementById('conv-stat')) {
                document.getElementById('conv-stat').textContent = `Displacement Computed: Solved in ${(performance.now() - t0).toFixed(1)}ms`;
            }
            activeState = JSON.parse(JSON.stringify(targetState));
        }
 catch (err) {
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Error: ${err.message}`;
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
