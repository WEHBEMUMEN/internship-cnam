/**
 * Phase 5.3a - UI and Simulation Manager
 * Coordinates sliders, ECSW toggle switches, benchmarking metrics, and Chart.js real-time displays.
 */

document.addEventListener('DOMContentLoaded', () => {
    // 1. Initialize Engines
    const engine = new window.NURBS2D();
    const solver = new window.MappingSolver(engine);
    
    // Initialize solvers reference state
    solver.initializeReference(0, 0);
    solver.mappedSolverType = document.getElementById('select-solver-type').value;

    // 2. Instantiate Visualizers
    const naiveVisuals = new window.VisualsEngine('naive-viewport', false);
    const mappedVisuals = new window.VisualsEngine('mapped-viewport', true);

    // Initial visuals update
    updateVisualizers();

    // 3. UI Elements References
    const sliderMu = document.getElementById('param-mu');
    const displayMu = document.getElementById('val-mu');
    const sliderR = document.getElementById('param-r');
    const displayR = document.getElementById('val-r');
    const sliderE = document.getElementById('param-E');
    const displayE = document.getElementById('val-E');
    const sliderNu = document.getElementById('param-nu');
    const displayNu = document.getElementById('val-nu');
    const sliderAmp = document.getElementById('param-amp');
    const displayAmp = document.getElementById('val-amp');
    const sliderFreq = document.getElementById('param-freq');
    const displayFreq = document.getElementById('val-freq');

    // Toggles
    const toggleControlNet = document.getElementById('toggle-cp');
    const toggleBCs = document.getElementById('toggle-bc');
    const toggleWarp = document.getElementById('toggle-warp');
    const selectViewMode = document.getElementById('select-viewmode');
    
    // ECSW Toggles & Displays
    const selectSolverType = document.getElementById('select-solver-type');
    const activeElementsToggleContainer = document.getElementById('active-elements-toggle-container');
    const toggleActiveElements = document.getElementById('toggle-active-elements');
    const ecswStatsContainer = document.getElementById('ecsw-stats-container');
    const ecswActiveCount = document.getElementById('ecsw-active-count');
    const ecswSparsity = document.getElementById('ecsw-sparsity');
    const mappedViewportLbl = document.getElementById('mapped-viewport-lbl');
    const statMappedLbl = document.getElementById('stat-mapped-lbl');

    // Controls
    const btnPlay = document.getElementById('btn-play');
    const btnReset = document.getElementById('btn-reset');
    const btnStep = document.getElementById('btn-step');

    // Stats Displays
    const displayDofs = document.getElementById('stat-dofs');
    const displayNaiveTime = document.getElementById('stat-naive-time');
    const displayMappedTime = document.getElementById('stat-mapped-time');
    const displaySpeedup = document.getElementById('stat-speedup');
    const displayFps = document.getElementById('stat-fps');

    // SVD Status Badge
    const svdNaiveBadge = document.getElementById('svd-naive-status');
    const svdMappedBadge = document.getElementById('svd-mapped-status');

    // 4. Set up Chart.js for real-time tip displacement comparison
    const ctx = document.getElementById('chart-comparison').getContext('2d');
    const comparisonChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [
                {
                    label: 'Naive Physical FOM (mm)',
                    data: [],
                    borderColor: '#06b6d4',
                    borderWidth: 2,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.3
                },
                {
                    label: 'Mapped Pullback FOM (mm)',
                    data: [],
                    borderColor: '#a78bfa',
                    borderWidth: 2,
                    borderDash: [5, 5],
                    pointRadius: 0,
                    fill: false,
                    tension: 0.3
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                x: {
                    grid: { color: 'rgba(255, 255, 255, 0.05)' },
                    ticks: { color: '#94a3b8', font: { family: 'Outfit', size: 10 } },
                    title: { display: true, text: 'Simulation Time (s)', color: '#94a3b8' }
                },
                y: {
                    grid: { color: 'rgba(255, 255, 255, 0.05)' },
                    ticks: { color: '#94a3b8', font: { family: 'Outfit', size: 10 } },
                    title: { display: true, text: 'Tip Vertical Displacement (mm)', color: '#94a3b8' }
                }
            },
            plugins: {
                legend: {
                    labels: {
                        color: '#f1f5f9',
                        font: { family: 'Outfit', size: 11 }
                    }
                }
            }
        }
    });

    // 5. State Management
    let isPlaying = false;
    let lastFrameTime = performance.now();
    let frameCount = 0;
    let fpsTimer = 0;
    let chartTimeWindow = 4.0; // Show last 4 seconds

    // Set DOFs stat
    const nU = solver.referencePatch.controlPoints.length;
    const nV = solver.referencePatch.controlPoints[0].length;
    displayDofs.innerText = `${nU * nV * 2} DOFs`;

    // 6. Register Event Listeners
    function registerSlider(slider, display, key, isSolver = true) {
        if (!slider) return;
        slider.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            display.innerText = val.toFixed(slider.step.includes('.') ? slider.step.split('.')[1].length : 0);
            
            if (isSolver) {
                solver[key] = val;
            }

            if (key === 'mu' || key === 'r') {
                updateSVDStatus();
            }

            if (!isPlaying) {
                updateVisualizers();
            }
        });
    }

    registerSlider(sliderMu, displayMu, 'mu');
    registerSlider(sliderR, displayR, 'r');
    registerSlider(sliderE, displayE, 'E');
    registerSlider(sliderNu, displayNu, 'nu');
    registerSlider(sliderAmp, displayAmp, 'forceAmp');
    registerSlider(sliderFreq, displayFreq, 'forceFreq');

    // ECSW Toggle
    // Solver Selection Dropdown
    if (selectSolverType) {
        selectSolverType.addEventListener('change', (e) => {
            solver.mappedSolverType = e.target.value;
            solver.useECSW = (solver.mappedSolverType === 'ECSW');
            
            if (solver.mappedSolverType === 'ECSW') {
                ecswStatsContainer.style.display = 'block';
                if (activeElementsToggleContainer) activeElementsToggleContainer.style.display = 'flex';
                ecswActiveCount.innerText = `${solver.ecswElements.length}/${solver.precomputedGauss.length}`;
                const ratio = ((solver.precomputedGauss.length - solver.ecswElements.length) / solver.precomputedGauss.length * 100).toFixed(1);
                ecswSparsity.innerText = `${ratio}%`;
                
                mappedViewportLbl.innerText = "ECSW Mapped Pullback Domain Ω̂";
                statMappedLbl.innerText = "ECSW Mapped Assembly";
                
                // Update chart labels
                comparisonChart.data.datasets[1].label = 'ECSW Mapped Pullback (mm)';
                comparisonChart.data.datasets[1].borderColor = '#10b981'; // Glowing green
            } else if (solver.mappedSolverType === 'Galerkin') {
                ecswStatsContainer.style.display = 'block';
                if (activeElementsToggleContainer) activeElementsToggleContainer.style.display = 'none';
                ecswActiveCount.innerText = `${solver.precomputedGauss.length}/${solver.precomputedGauss.length}`;
                ecswSparsity.innerText = '0.0%';
                
                mappedViewportLbl.innerText = "Galerkin Mapped Pullback Domain Ω̂";
                statMappedLbl.innerText = "Galerkin ROM Assembly";
                
                // Update chart labels
                comparisonChart.data.datasets[1].label = 'Galerkin Mapped Pullback (mm)';
                comparisonChart.data.datasets[1].borderColor = '#f59e0b'; // Orange
            } else {
                ecswStatsContainer.style.display = 'none';
                if (activeElementsToggleContainer) activeElementsToggleContainer.style.display = 'none';
                mappedViewportLbl.innerText = "Mapped Pullback Domain Ω̂";
                statMappedLbl.innerText = "Mapped Online Assembly";
                
                comparisonChart.data.datasets[1].label = 'Mapped Pullback FOM (mm)';
                comparisonChart.data.datasets[1].borderColor = '#a78bfa'; // Purple
            }
            comparisonChart.update('none');

            if (!isPlaying) {
                updateVisualizers();
            }
        });
    }

    // Toggle active elements highlighting
    if (toggleActiveElements) {
        toggleActiveElements.addEventListener('change', (e) => {
            mappedVisuals.showActiveElements = e.target.checked;
            if (!isPlaying) {
                updateVisualizers();
            }
        });
    }

    // Controls
    if (btnPlay) {
        btnPlay.addEventListener('click', () => {
            isPlaying = !isPlaying;
            if (isPlaying) {
                btnPlay.innerHTML = '<i class="fa-solid fa-pause"></i> Pause Sim';
                btnPlay.style.background = 'linear-gradient(135deg, #ef4444, #b91c1c)';
                btnPlay.style.boxShadow = '0 10px 20px -5px rgba(239, 68, 68, 0.4)';
                lastFrameTime = performance.now();
                animate();
            } else {
                btnPlay.innerHTML = '<i class="fa-solid fa-play"></i> Run Simulation';
                btnPlay.style.background = 'linear-gradient(135deg, var(--primary), var(--primary-dark))';
                btnPlay.style.boxShadow = '0 10px 20px -5px rgba(124, 58, 237, 0.4)';
            }
        });
    }

    if (btnReset) {
        btnReset.addEventListener('click', () => {
            solver.resetDynamicStates();
            comparisonChart.data.labels = [];
            comparisonChart.data.datasets[0].data = [];
            comparisonChart.data.datasets[1].data = [];
            comparisonChart.update();
            updateVisualizers();
            if (displayNaiveTime) displayNaiveTime.innerText = '--';
            if (displayMappedTime) displayMappedTime.innerText = '--';
            if (displaySpeedup) {
                displaySpeedup.innerText = '--';
                displaySpeedup.className = 'metric-card';
            }
        });
    }

    if (btnStep) {
        btnStep.addEventListener('click', () => {
            runSingleStep();
        });
    }

    // Toggle events
    if (toggleControlNet) {
        toggleControlNet.addEventListener('change', (e) => {
            naiveVisuals.showControlNet = e.target.checked;
            mappedVisuals.showControlNet = e.target.checked;
            updateVisualizers();
        });
    }

    if (toggleBCs) {
        toggleBCs.addEventListener('change', (e) => {
            naiveVisuals.showBCs = e.target.checked;
            mappedVisuals.showBCs = e.target.checked;
            updateVisualizers();
        });
    }

    if (toggleWarp) {
        toggleWarp.addEventListener('change', (e) => {
            mappedVisuals.warpToPhysical = e.target.checked;
            updateVisualizers();
        });
    }

    if (selectViewMode) {
        selectViewMode.addEventListener('change', (e) => {
            naiveVisuals.viewMode = e.target.value;
            mappedVisuals.viewMode = e.target.value;
            updateVisualizers();
        });
    }

    // Initialize solver type UI state
    if (selectSolverType) {
        selectSolverType.dispatchEvent(new Event('change'));
    }

    // Initialize SVD indicators
    updateSVDStatus();

    // 7. Core Loop & Solves
    function runSingleStep() {
        const result = solver.stepDynamics(solver.mu, solver.r);

        const tipNodeIdx = ((nU - 1) * 3 + 1) * 2 + 1;
        const dispYNaive = result.uNaive[tipNodeIdx];
        const dispYMapped = result.uMapped[tipNodeIdx];

        // Update charts
        const curTime = result.time;
        comparisonChart.data.labels.push(curTime.toFixed(2));
        comparisonChart.data.datasets[0].data.push(dispYNaive);
        comparisonChart.data.datasets[1].data.push(dispYMapped);

        if (comparisonChart.data.labels.length > Math.round(chartTimeWindow / solver.dt)) {
            comparisonChart.data.labels.shift();
            comparisonChart.data.datasets[0].data.shift();
            comparisonChart.data.datasets[1].data.shift();
        }
        comparisonChart.update('none');

        // Benchmarks display
        if (displayNaiveTime) displayNaiveTime.innerText = `${result.tNaiveAssembly.toFixed(2)} ms`;
        
        // Show assembly time in microseconds if really fast, or ms otherwise
        if (displayMappedTime) {
            if (result.tMappedAssembly < 1.0) {
                displayMappedTime.innerText = `${(result.tMappedAssembly * 1000).toFixed(0)} µs`;
            } else {
                displayMappedTime.innerText = `${result.tMappedAssembly.toFixed(2)} ms`;
            }
        }

        if (displaySpeedup) {
            const speedupVal = result.tNaiveAssembly / result.tMappedAssembly;
            displaySpeedup.innerText = `${speedupVal.toFixed(1)}x`;
            
            if (speedupVal > 30) {
                displaySpeedup.className = 'stat-val speedup-extreme';
            } else if (speedupVal > 10) {
                displaySpeedup.className = 'stat-val speedup-good';
            } else {
                displaySpeedup.className = 'stat-val speedup-low';
            }
        }

        updateVisualizers();
    }

    function updateVisualizers() {
        const naivePatch = window.GeometryFactory.refine(
            window.GeometryFactory.generateNotchedBeam(solver.L, solver.H, solver.mu, solver.r),
            0, 0
        );

        naiveVisuals.update(solver, naivePatch, solver.uNaive);
        mappedVisuals.update(solver, solver.referencePatch, solver.uMapped, naivePatch);
    }

    function updateSVDStatus() {
        if (svdNaiveBadge) {
            svdNaiveBadge.className = 'svd-badge incompatible';
            svdNaiveBadge.innerHTML = '<i class="fa-solid fa-triangle-exclamation"></i> INCOMPATIBLE (DOFs vary)';
        }

        if (svdMappedBadge) {
            svdMappedBadge.className = 'svd-badge compatible';
            svdMappedBadge.innerHTML = '<i class="fa-solid fa-circle-check"></i> SVD COMPATIBLE (Static Grid)';
        }
    }

    function animate() {
        if (!isPlaying) return;

        runSingleStep();

        const now = performance.now();
        frameCount++;
        fpsTimer += (now - lastFrameTime);
        lastFrameTime = now;

        if (fpsTimer >= 1000) {
            if (displayFps) displayFps.innerText = Math.round((frameCount * 1000) / fpsTimer);
            frameCount = 0;
            fpsTimer = 0;
        }

        requestAnimationFrame(animate);
    }
});
