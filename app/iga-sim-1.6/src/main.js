import { NURBSEngine } from './engine/iga.js';
import { PhysicsEngine } from './engine/physics.js';
import { BasisPlot } from './ui/plot.js';

class MechanicsApp {
    constructor() {
        this.canvas = document.getElementById('plot-canvas');
        this.ctx = this.canvas.getContext('2d');
        
        // State
        this.degree = 3;
        this.numElements = 4;
        this.loadPos = 0.5;
        this.loadMag = 10.0;
        this.visualScale = 100.0;
        this.isStatics = true;
        this.isNonlinear = false;
        this.showFEM = false;
        this.bcLeft = true;
        this.bcRight = true;
        
        this.time = 0;
        this.animSpeed = 1.0;
        this.deflection = null;
        this.femDeflection = null;
        this.femDeflection = null;
        this.femQuality = 9;
        this.femType = 'linear';
        this.showAnalytical = false;
        
        // Initialize Engines
        this.initEngines();
        this.plot = new BasisPlot('plot-canvas');
        
        this.setupEventListeners();
        this.resize();
        window.addEventListener('resize', () => this.resize());
        
        this.updatePhysics();
        this.renderMath();
        this.animate();
    }

    renderMath() {
        const formula = document.getElementById('math-formula');
        if (!formula || !window.katex) return;
        
        let tex = this.isNonlinear 
            ? "[K(u)]\\{u\\} = \\{F(u)\\}" 
            : "[K]\\{u\\} = \\{F\\}";
        
        window.katex.render(tex, formula, {
            throwOnError: false,
            displayMode: true
        });
    }

    initEngines() {
        // Create a simple straight line NURBS as our structure
        const numCP = this.numElements + this.degree;
        const knots = [];
        for (let i = 0; i <= this.degree; i++) knots.push(0);
        for (let i = 1; i < this.numElements; i++) knots.push(i / this.numElements);
        for (let i = 0; i <= this.degree; i++) knots.push(1);

        const controlPoints = [];
        for(let i = 0; i < numCP; i++) {
            controlPoints.push({ x: i / (numCP-1), y: 0.5, w: 1.0 });
        }

        this.nurbs = new NURBSEngine(this.degree, knots, controlPoints);
        this.physics = new PhysicsEngine(this.nurbs);
    }

    updatePhysics() {
        const numCP = this.nurbs.controlPoints.length;
        
        // Exact Point Load Assembly
        this.loadF = this.physics.assembleIGALoad(this.loadPos, this.loadMag);

        // IGA BCs: Clamped means pinning first/last TWO control points (locks value AND derivative)
        const igaBCs = [];
        if (this.bcLeft) {
            igaBCs.push({ index: 0, value: 0 });
            igaBCs.push({ index: 1, value: 0 });
        }
        if (this.bcRight) {
            igaBCs.push({ index: numCP - 1, value: 0 });
            igaBCs.push({ index: numCP - 2, value: 0 });
        }

        // Always solve the static shape — in dynamics mode it serves as the vibration mode shape
        if (this.isNonlinear) {
            this.deflection = this.physics.solveNonlinear(this.loadF, igaBCs);
        } else {
            this.deflection = this.physics.solveStatics(this.loadF, igaBCs);
        }

        // FEM BCs: Hermite elements use DOF indices 0 and last node
        const femBCs = [];
        if (this.bcLeft) femBCs.push({ index: 0, value: 0 });
        if (this.bcRight) femBCs.push({ index: numCP - 1, value: 0 });

        if (this.showFEM) {
            if (this.isNonlinear) {
                this.femDeflection = this.physics.solveNonlinearFEM(this.loadPos, this.loadMag, femBCs, this.femQuality, this.femType);
            } else {
                this.femDeflection = this.physics.solveFEM(this.loadPos, this.loadMag, femBCs, this.femQuality, this.femType);
            }
        }
    }

    setupEventListeners() {
        const inputs = {
            'input-degree': (v) => { this.degree = parseInt(v); document.getElementById('degree-val').textContent = v; this.initEngines(); },
            'input-elements': (v) => { this.numElements = parseInt(v); document.getElementById('elements-val').textContent = v; this.initEngines(); },
            'input-load-pos': (v) => { this.loadPos = parseFloat(v); document.getElementById('load-pos-val').textContent = v; },
            'input-load-mag': (v) => { this.loadMag = parseFloat(v); document.getElementById('load-mag-val').textContent = v; },
            'input-fem-quality': (v) => { this.femQuality = parseInt(v); document.getElementById('fem-quality-val').textContent = v; }
        };

        Object.entries(inputs).forEach(([id, fn]) => {
            document.getElementById(id).addEventListener('input', (e) => {
                fn(e.target.value);
                this.updatePhysics();
            });
        });

        document.getElementById('toggle-nonlinear').addEventListener('change', (e) => {
            this.isNonlinear = e.target.checked;
            this.updatePhysics();
            this.renderMath();
        });

        // Moved to the block below to handle legend logic

        document.getElementById('bc-left').addEventListener('click', (e) => {
            this.bcLeft = !this.bcLeft;
            e.target.classList.toggle('active', this.bcLeft);
            e.target.textContent = `Left: ${this.bcLeft ? 'Fixed' : 'Free'}`;
            this.updatePhysics();
        });

        document.getElementById('bc-right').addEventListener('click', (e) => {
            this.bcRight = !this.bcRight;
            e.target.classList.toggle('active', this.bcRight);
            e.target.textContent = `Right: ${this.bcRight ? 'Fixed' : 'Free'}`;
            this.updatePhysics();
        });

        document.getElementById('mode-statics').addEventListener('click', (e) => {
            this.isStatics = true;
            e.target.classList.add('active');
            document.getElementById('mode-dynamics').classList.remove('active');
            document.getElementById('dynamics-settings').style.display = 'none';
            document.getElementById('state-desc').textContent = 'Solving natively on NURBS basis for C¹-smooth continuous deflection.';
            this.updatePhysics();
            this.renderMath();
        });

        document.getElementById('mode-dynamics').addEventListener('click', (e) => {
            this.isStatics = false;
            this.time = 0;
            e.target.classList.add('active');
            document.getElementById('mode-statics').classList.remove('active');
            document.getElementById('dynamics-settings').style.display = 'block';
            document.getElementById('state-desc').textContent = 'Simulating structural vibration mode shapes using NURBS Basis.';
            this.updatePhysics();
            this.renderMath();
        });

        document.getElementById('input-speed').addEventListener('input', (e) => {
            this.animSpeed = parseFloat(e.target.value);
            document.getElementById('speed-val').textContent = this.animSpeed.toFixed(1);
        });

        document.getElementById('btn-reset').addEventListener('click', () => location.reload());
        
        // View modes
        document.getElementById('view-basis').addEventListener('click', (e) => {
            e.target.classList.add('active');
            document.getElementById('view-deflection').classList.remove('active');
            document.getElementById('canvas-legend').style.display = 'none';
        });
        document.getElementById('view-deflection').addEventListener('click', (e) => {
            e.target.classList.add('active');
            document.getElementById('view-basis').classList.remove('active');
            document.getElementById('canvas-legend').style.display = 'flex';
        });

        // Legend logic for FEM
        document.getElementById('toggle-fem').addEventListener('change', (e) => {
            this.showFEM = e.target.checked;
            document.getElementById('fem-settings-container').style.display = this.showFEM ? 'block' : 'none';
            document.getElementById('legend-fem').style.display = this.showFEM ? 'flex' : 'none';
            this.updatePhysics();
        });

        // Legend logic for Analytical Beam Theory
        document.getElementById('toggle-analytical').addEventListener('change', (e) => {
            this.showAnalytical = e.target.checked;
            document.getElementById('legend-analytical').style.display = this.showAnalytical ? 'flex' : 'none';
            this.updatePhysics();
        });

        document.getElementById('input-scale-factor').addEventListener('input', (e) => {
            this.visualScale = parseFloat(e.target.value);
            document.getElementById('scale-factor-val').textContent = this.visualScale.toString();
        });

        // Benchmark Modal Logic
        document.getElementById('btn-benchmark').addEventListener('click', () => {
            const results = this.physics.runBenchmark(this.loadMag, 100);
            const report = `===========================================
    PERFORMANCE & CONVERGENCE BENCHMARK
===========================================

System Information:
-------------------
IGA Degree (p) : ${this.degree}
IGA Elements   : ${this.numElements}
IGA Total DOFs : ${results.dofs}

FEM Elements   : ${Math.floor(results.femDofs/2)} (Quadratic) / ${results.femDofs-1} (Linear)
FEM Total DOFs : ${results.femDofs}

Average Solve Times (100 Iterations):
-------------------------------------
[1] Isogeometric Analysis Solver
    Time: ${results.iga} ms

[2] Linear FEM Solver (2-node)
    Time: ${results.lin} ms

[3] Quadratic FEM Solver (3-node)
    Time: ${results.quad} ms
    
Conclusion:
-----------
IGA enables superior smoothness and accuracy 
using fewer degrees of freedom compared to 
traditional Classical FEM approaches.`;
            
            document.getElementById('benchmark-results').textContent = report;
            document.getElementById('modal-container').style.display = 'flex';
        });

        document.getElementById('btn-close-modal').addEventListener('click', () => {
            document.getElementById('modal-container').style.display = 'none';
        });

        // Toolbar
        document.getElementById('btn-zoom-in').addEventListener('click', () => this.plot.handleZoom(1.5));
        document.getElementById('btn-zoom-out').addEventListener('click', () => this.plot.handleZoom(0.7));
        document.getElementById('btn-zoom-reset').addEventListener('click', () => this.plot.resetView());
        
        // Default legend setup
        document.getElementById('canvas-legend').style.display = 'flex';
    }

    resize() {
        this.canvas.width = this.canvas.parentElement.clientWidth;
        this.canvas.height = this.canvas.parentElement.clientHeight;
    }

    drawLoadVector() {
        const pt = this.plot.worldToScreen(this.loadPos, 0.5);
        const length = this.loadMag * this.plot.camera.zoom * 1.5;
        
        this.ctx.strokeStyle = '#ef4444';
        this.ctx.lineWidth = 3;
        this.ctx.beginPath();
        this.ctx.moveTo(pt.x, pt.y);
        this.ctx.lineTo(pt.x, pt.y + length);
        this.ctx.stroke();
        
        // Arrowhead
        const head = length > 0 ? 10 : -10;
        this.ctx.beginPath();
        this.ctx.moveTo(pt.x - 5, pt.y + length - head);
        this.ctx.lineTo(pt.x, pt.y + length);
        this.ctx.lineTo(pt.x + 5, pt.y + length - head);
        this.ctx.stroke();
    }

    animate() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.plot.drawGrid();
        
        const isBasisView = document.getElementById('view-basis').classList.contains('active');
        
        if (isBasisView) {
            this.plot.drawBasis(this.nurbs);
        } else {
            // Draw Deflection Shape
            let dynamicScale = 1.0;
            if (!this.isStatics) {
                this.time += 0.05 * this.animSpeed;
                dynamicScale = Math.sin(this.time);
            }

            const renderedScale = this.visualScale * 0.005 * dynamicScale;

            // Create temporary NURBS with Y coordinates shifted by deflection
            // Sign mapping: Y increases upwards in math, so downwards physical deflection means we subtract
            const shiftedCPs = this.nurbs.controlPoints.map((cp, i) => ({
                x: cp.x,
                y: 0.5 - (this.deflection[i] * renderedScale), 
                w: cp.w
            }));
            
            const tempNurbs = new NURBSEngine(this.nurbs.degree, this.nurbs.knots, shiftedCPs);
            this.plot.drawCurve(tempNurbs);
            
            // Draw FEM Comparison if enabled
            if (this.showFEM && this.femDeflection) {
                this.plot.drawFEMCurve(this.femDeflection, renderedScale);
            }

            // Draw Analytical Theory if enabled
            if (this.showAnalytical) {
                this.plot.drawAnalyticalCurve({
                    loadPos: this.loadPos,
                    effectiveLoad: this.loadMag,
                    engine: this.physics,
                    visualScale: renderedScale
                });
            }

            this.drawLoadVector();
            
            // Handle hover logic
            const mX = this.plot.mousePos.worldX;
            if (mX >= 0 && mX <= 1) {
                // Approximate exact value from NURBS or array
                const numCP = this.deflection.length;
                const idx = Math.max(0, Math.min(numCP - 1, Math.round(mX * (numCP - 1))));
                const val = this.deflection[idx];
                this.plot.hoverValue = val;
            } else {
                this.plot.hoverValue = null;
            }
        }

        this.plot.drawCrosshair();
        requestAnimationFrame(() => this.animate());
    }
}

new MechanicsApp();
