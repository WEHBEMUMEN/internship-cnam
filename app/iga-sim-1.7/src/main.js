import { NURBSEngine } from './engine/iga.js';
import { PhysicsEngine, ReferenceFEM } from './engine/physics.js';
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
        this.showFEM = true; // Enabled by default for 1.7 Bench
        this.bcLeft = true;
        this.bcRight = true;
        
        this.time = 0;
        this.animSpeed = 1.0;
        this.isDamping = true;
        this.dynamicAmplitude = 1.0;
        this.deflection = null;
        this.femDeflection = null;
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
        
        // For Debugging
        window.app = this;
    }

    renderMath() {
        const formula = document.getElementById('math-formula');
        if (!formula) return;
        
        let tex = this.isNonlinear 
            ? "[K(u)]\\{u\\} = \\{F(u)\\}" 
            : "[K]\\{u\\} = \\{F\\}";
        
        if (window.katex) {
            window.katex.render(tex, formula, {
                throwOnError: false,
                displayMode: true
            });
        } else {
            formula.textContent = tex;
        }
    }

    initEngines() {
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
        this.referenceFEM = new ReferenceFEM(100);
    }

    updatePhysics() {
        const numCP = this.nurbs.controlPoints.length;
        this.loadF = this.physics.assembleIGALoad(this.loadPos, this.loadMag);

        const igaBCs = [];
        if (this.bcLeft) {
            igaBCs.push({ index: 0, value: 0 });
            igaBCs.push({ index: 1, value: 0 });
        }
        if (this.bcRight) {
            igaBCs.push({ index: numCP - 1, value: 0 });
            igaBCs.push({ index: numCP - 2, value: 0 });
        }

        this.deflection = this.physics.solveStatics(this.loadF, igaBCs);
        
        // Unified Nonlinear Softening (Applied to both ROM and Reference for parity)
        const beta = 25.0;
        if (this.isNonlinear) {
            this.deflection = this.deflection.map(val => val * (1.0 + beta * Math.pow(val, 2)));
        }

        if (this.showFEM) {
            // Apply similar logic for Reference FEM (Clamping)
            const femBCs = [];
            if (this.bcLeft) {
                femBCs.push({ index: 0, value: 0 });
                femBCs.push({ index: 1, value: 0 });
            }
            if (this.bcRight) {
                const lastNodeDof = 202 - 2; // (100+1)*2 = 202. last node starts at 200
                femBCs.push({ index: lastNodeDof, value: 0 });
                femBCs.push({ index: lastNodeDof + 1, value: 0 });
            }
            this.femDeflection = this.referenceFEM.solve(this.loadPos, this.loadMag, femBCs);
            
            if (this.isNonlinear) {
                this.femDeflection = this.femDeflection.map(p => ({
                    x: p.x,
                    y: p.y * (1.0 + beta * Math.pow(p.y, 2))
                }));
            }
        }
    }

    setupEventListeners() {
        const inputMap = {
            'input-degree': (v) => { this.degree = parseInt(v); document.getElementById('degree-val').textContent = v; this.initEngines(); },
            'input-elements': (v) => { this.numElements = parseInt(v); document.getElementById('elements-val').textContent = v; this.initEngines(); },
            'input-load-pos': (v) => { this.loadPos = parseFloat(v); document.getElementById('load-pos-val').textContent = v; },
            'input-load-mag': (v) => { this.loadMag = parseFloat(v); document.getElementById('load-mag-val').textContent = v; },
            'input-scale-factor': (v) => { this.visualScale = parseFloat(v); document.getElementById('scale-factor-val').textContent = v; }
        };

        Object.entries(inputMap).forEach(([id, fn]) => {
            const el = document.getElementById(id);
            if (el) el.addEventListener('input', (e) => {
                fn(e.target.value);
                this.updatePhysics();
            });
        });

        const toggles = {
            'toggle-nonlinear': (v) => { this.isNonlinear = v; this.updatePhysics(); this.renderMath(); },
            'toggle-fem': (v) => { this.showFEM = v; document.getElementById('fem-settings-container').style.display = v ? 'block' : 'none'; document.getElementById('legend-fem').style.display = v ? 'flex' : 'none'; this.updatePhysics(); },
            'toggle-analytical': (v) => { this.showAnalytical = v; document.getElementById('legend-analytical').style.display = v ? 'flex' : 'none'; this.updatePhysics(); },
            'toggle-damping': (v) => { this.isDamping = v; }
        };

        Object.entries(toggles).forEach(([id, fn]) => {
            const el = document.getElementById(id);
            if (el) el.addEventListener('change', (e) => fn(e.target.checked));
        });

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
            this.updatePhysics();
        });

        document.getElementById('mode-dynamics').addEventListener('click', (e) => {
            this.isStatics = false;
            this.time = 0;
            this.dynamicAmplitude = 1.0;
            e.target.classList.add('active');
            document.getElementById('mode-statics').classList.remove('active');
            document.getElementById('dynamics-settings').style.display = 'block';
            this.updatePhysics();
        });

        document.getElementById('btn-reset').addEventListener('click', () => location.reload());
        
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

        document.getElementById('btn-benchmark').addEventListener('click', () => {
            const results = this.physics.runBenchmark(this.loadMag, 100);
            const report = `PERFORMANCE BENCHMARK\nIGA DOFs: ${results.dofs}\nSolve Time: ${results.iga} ms\nFEM (Ref) DOFs: 202`;
            document.getElementById('benchmark-results').textContent = report;
            document.getElementById('modal-container').style.display = 'flex';
        });

        document.getElementById('btn-close-modal').addEventListener('click', () => {
            document.getElementById('modal-container').style.display = 'none';
        });
    }

    resize() {
        if (!this.canvas) return;
        this.canvas.width = this.canvas.parentElement.clientWidth;
        this.canvas.height = this.canvas.parentElement.clientHeight;
    }

    drawLoadVector(renderedScale) {
        const pt = this.plot.worldToScreen(this.loadPos, 0.5);
        const length = this.loadMag * this.plot.camera.zoom * 1.5;
        this.ctx.strokeStyle = '#ef4444';
        this.ctx.lineWidth = 3;
        this.ctx.beginPath();
        this.ctx.moveTo(pt.x, pt.y);
        this.ctx.lineTo(pt.x, pt.y + length);
        this.ctx.stroke();
    }

    animate() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.plot.drawGrid();
        
        const isBasisView = document.getElementById('view-basis').classList.contains('active');
        
        if (isBasisView) {
            this.plot.drawBasis(this.nurbs);
        } else {
            let dynamicScale = 1.0;
            if (!this.isStatics) {
                this.time += 0.05 * this.animSpeed;
                if (this.isDamping) {
                    this.dynamicAmplitude *= 0.99;
                } else {
                    this.dynamicAmplitude = 1.0;
                }
                dynamicScale = Math.sin(this.time) * this.dynamicAmplitude;
            }

            const renderedScale = this.visualScale * 0.005 * dynamicScale;

            if (this.deflection) {
                const shiftedCPs = this.nurbs.controlPoints.map((cp, i) => ({
                    x: cp.x,
                    y: 0.5 - (this.deflection[i] * renderedScale), 
                    w: cp.w
                }));
                const tempNurbs = new NURBSEngine(this.nurbs.degree, this.nurbs.knots, shiftedCPs);
                this.plot.drawCurve(tempNurbs, '#3b82f6');
                
                if (this.showFEM && this.femDeflection) {
                    this.plot.drawReferenceFEM(this.femDeflection, renderedScale);
                }
                
                if (this.showAnalytical) {
                    this.plot.drawAnalyticalCurve({
                        loadPos: this.loadPos,
                        effectiveLoad: this.loadMag,
                        engine: this.physics,
                        visualScale: renderedScale
                    });
                }
                this.updateROMStats();
            }
            this.drawLoadVector(renderedScale);
        }

        this.plot.drawCrosshair();
        requestAnimationFrame(() => this.animate());
    }

    updateROMStats() {
        const romDofs = this.physics.getDegreesOfFreedom();
        const refDofs = this.referenceFEM ? this.referenceFEM.getDegreesOfFreedom() : 202;
        const savings = ((1 - romDofs / refDofs) * 100).toFixed(1);

        const romDofEl = document.getElementById('rom-dofs');
        const refDofEl = document.getElementById('ref-dofs');
        const savingsEl = document.getElementById('rom-savings');

        if (romDofEl) romDofEl.textContent = romDofs;
        if (refDofEl) refDofEl.textContent = refDofs;
        if (savingsEl) savingsEl.textContent = savings + '%';
    }
}

new MechanicsApp();
