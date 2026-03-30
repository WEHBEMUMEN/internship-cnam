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
        this.stiffness = 100.0;
        this.bcLeft = true;
        this.bcRight = true;
        
        this.deflection = null;
        this.romDeflection = null;
        this.torqueResult = null;
        this.isTorqueMode = false;
        this.isNonlinear = true; // Exclusively Non-Linear
        this.physicsMode = 'bending';
        
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
        
        let tex = this.isTorqueMode
            ? "\\kappa = M / EI"
            : "[K(u)]\\{u\\} = \\{F(u)\\}";
        
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
        this.physics.E = this.stiffness;
        this.referenceFEM = new ReferenceFEM(40);
        this.referenceFEM.E = this.stiffness;
    }

    updatePhysics() {
        if (this.isTorqueMode) {
            this.torqueResult = this.physics.solveTorqueToCircle(this.loadMag);
            this.deflection = this.romDeflection = this.femDeflection = null;
            this.updateROMStats(true);
            this.renderMath();
            return;
        }
        this.torqueResult = null;

        const numCP = this.nurbs.controlPoints.length;
        this.physics.physicsMode = this.physicsMode;
        
        // Apply BCs
        const igaBCs = [];
        if (this.bcLeft) {
            igaBCs.push({ index: 0, value: 0 });
            if (this.physicsMode === 'bending') igaBCs.push({ index: 1, value: 0 });
        }
        if (this.bcRight) {
            igaBCs.push({ index: numCP - 1, value: 0 });
            if (this.physicsMode === 'bending') igaBCs.push({ index: numCP - 2, value: 0 });
        }

        const loadF = this.physics.assembleIGALoad(this.loadPos, this.loadMag);
        
        // Non-Linear Solver Selection
        const numModes = Math.min(3, numCP - 2);
        this.romModes = numModes;

        this.deflection = this.physics.solveNonLinearStatics(loadF, igaBCs, 10);
        this.romDeflection = this.physics.solveNonLinearROM(loadF, igaBCs, this.romModes, 10);

        this.updateROMStats();
        this.renderMath();
    }

    setupEventListeners() {
        const inputMap = {
            'input-degree': (v) => { this.degree = parseInt(v); document.getElementById('degree-val').textContent = v; this.initEngines(); },
            'input-elements': (v) => { this.numElements = parseInt(v); document.getElementById('elements-val').textContent = v; this.initEngines(); },
            'input-stiffness': (v) => { 
                this.stiffness = parseFloat(v); 
                document.getElementById('stiffness-val').textContent = v;
                this.physics.E = this.stiffness;
            },
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

        document.getElementById('select-physics-mode').addEventListener('change', (e) => {
            this.physicsMode = e.target.value;
            if (this.physicsMode === 'axial' && this.bcRight) {
                this.bcRight = false;
                const bcBtn = document.getElementById('bc-right');
                bcBtn.classList.remove('active');
                bcBtn.textContent = 'Right: Free';
            }
            this.updatePhysics();
        });

        document.getElementById('btn-torque-circle').addEventListener('click', () => {
            this.isTorqueMode = !this.isTorqueMode;
            document.getElementById('btn-torque-circle').classList.toggle('active', this.isTorqueMode);
            
            const torqueModeActive = this.isTorqueMode;
            document.getElementById('select-physics-mode').disabled = torqueModeActive;
            document.getElementById('input-load-pos').disabled = torqueModeActive;
            document.getElementById('bc-left').disabled = torqueModeActive;
            document.getElementById('bc-right').disabled = torqueModeActive;

            this.updatePhysics();
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
    }

    resize() {
        if (!this.canvas) return;
        this.canvas.width = this.canvas.parentElement.clientWidth;
        this.canvas.height = this.canvas.parentElement.clientHeight;
    }

    drawLoadVector(renderedScale) {
        if (this.isTorqueMode) {
            const p0 = this.plot.worldToScreen(0, 0.5);
            const p1 = this.plot.worldToScreen(1, 0.5);
            this.ctx.strokeStyle = '#a78bfa';
            this.ctx.lineWidth = 2;
            [p0, p1].forEach(p => {
                this.ctx.beginPath();
                this.ctx.arc(p.x, p.y, 15, -Math.PI * 0.8, Math.PI * 0.8, this.loadMag < 0);
                this.ctx.stroke();
            });
            return;
        }

        const currentLoadPos = this.physicsMode === 'axial' ? 1.0 : this.loadPos;
        const pt = this.plot.worldToScreen(currentLoadPos, 0.5);
        const length = (this.loadMag / 5.0) * this.plot.camera.zoom; 
        
        this.ctx.strokeStyle = '#ef4444';
        this.ctx.lineWidth = 3;
        this.ctx.beginPath();
        
        this.ctx.moveTo(pt.x, pt.y);
        if (this.physicsMode === 'axial') {
            this.ctx.lineTo(pt.x + length, pt.y);
            const dir = length > 0 ? 1 : -1;
            this.ctx.moveTo(pt.x + length, pt.y);
            this.ctx.lineTo(pt.x + length - (10 * dir), pt.y - 5);
            this.ctx.moveTo(pt.x + length, pt.y);
            this.ctx.lineTo(pt.x + length - (10 * dir), pt.y + 5);
        } else {
            this.ctx.lineTo(pt.x, pt.y + length);
            const dir = length > 0 ? 1 : -1;
            this.ctx.moveTo(pt.x, pt.y + length);
            this.ctx.lineTo(pt.x - 5, pt.y + length - (10 * dir));
            this.ctx.moveTo(pt.x, pt.y + length);
            this.ctx.lineTo(pt.x + 5, pt.y + length - (10 * dir));
        }
        this.ctx.stroke();
    }

    animate() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.plot.drawGrid();
        
        const isBasisView = document.getElementById('view-basis').classList.contains('active');
        
        if (isBasisView) {
            this.plot.drawBasis(this.nurbs);
        } else if (this.isTorqueMode && this.torqueResult && this.torqueResult.length > 0) {
            this.ctx.strokeStyle = '#a78bfa';
            this.ctx.lineWidth = 3;
            this.ctx.beginPath();
            this.torqueResult.forEach((p, i) => {
                const screen = this.plot.worldToScreen(p.x, p.y);
                if (i === 0) this.ctx.moveTo(screen.x, screen.y);
                else this.ctx.lineTo(screen.x, screen.y);
            });
            this.ctx.stroke();
            this.drawLoadVector();
        } else {
            const renderedScale = this.visualScale * 0.005;

            if (this.deflection) {
                const shiftedROM = this.nurbs.controlPoints.map((cp, i) => {
                    const disp = this.romDeflection ? this.romDeflection[i] : 0;
                    if (this.physicsMode === 'axial') return { x: cp.x + (disp * renderedScale), y: 0.5, w: cp.w };
                    return { x: cp.x, y: 0.5 - (disp * renderedScale), w: cp.w };
                });
                const romCurve = new NURBSEngine(this.nurbs.degree, this.nurbs.knots, shiftedROM);
                this.plot.drawCurve(romCurve, 'rgba(59, 130, 246, 0.4)');

                const shiftedCPs = this.nurbs.controlPoints.map((cp, i) => {
                    if (this.physicsMode === 'axial') return { x: cp.x + (this.deflection[i] * renderedScale), y: 0.5, w: cp.w };
                    return { x: cp.x, y: 0.5 - (this.deflection[i] * renderedScale), w: cp.w };
                });
                const tempNurbs = new NURBSEngine(this.nurbs.degree, this.nurbs.knots, shiftedCPs);
                this.plot.drawCurve(tempNurbs, '#f8fafc');
            }

            this.drawLoadVector(renderedScale);
        }

        this.plot.drawCrosshair();
        requestAnimationFrame(() => this.animate());
    }

    updateROMStats(isTorqueMode = false) {
        const romDofEl = document.getElementById('rom-dofs');
        const refDofEl = document.getElementById('ref-dofs');
        const savingsEl = document.getElementById('rom-savings');
        const stateDesc = document.getElementById('state-desc');

        if (isTorqueMode) {
            if (romDofEl) romDofEl.textContent = 'N/A';
            if (refDofEl) refDofEl.textContent = 'N/A';
            if (savingsEl) savingsEl.textContent = 'N/A';
            if (stateDesc) stateDesc.textContent = 'Geometric non-linearity: moment applied to beam ends.';
            return;
        }

        const fullDofs = this.physics.getDegreesOfFreedom();
        const romDofs = this.romModes || 3;
        const savings = fullDofs > 0 ? ((1 - romDofs / fullDofs) * 100).toFixed(1) : '0.0';

        if (romDofEl) romDofEl.textContent = romDofs;
        if (refDofEl) refDofEl.textContent = fullDofs;
        if (savingsEl) savingsEl.textContent = savings + '%';
        
        if (stateDesc) stateDesc.textContent = `Reducing Full NL-IGA (${fullDofs} CP) to a ${romDofs}-mode NL-ROM. Gain: ${savings}%.`;
    }
}

new MechanicsApp();
