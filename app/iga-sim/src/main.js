import { KnotVector, BasisFunctions, NURBS, Curve } from './engine/iga.js';
import { BasisPlot } from './ui/plot.js';

class App {
    constructor() {
        this.p = 2; // Degree
        this.n = 6; // Number of control points
        this.mode = 'basis'; // 'basis' or 'curve'
        this.knotVector = new KnotVector(this.n, this.p);
        this.plot = new BasisPlot('sim-canvas');
        
        // NURBS State
        this.points = this.initPoints();
        this.weights = new Array(this.n).fill(1.0);
        
        // Interaction State
        this.draggingPoint = -1;
        
        this.setupEventListeners();
        this.updateKnotUI();
        this.updateWeightsUI();
        this.render();

        window.addEventListener('resize', () => this.render());
    }

    initPoints() {
        const pts = [];
        for (let i = 0; i < this.n; i++) {
            pts.push({
                x: i / (this.n - 1),
                y: 0.5 + 0.2 * Math.sin(i * Math.PI / 2)
            });
        }
        return pts;
    }

    setupEventListeners() {
        const degreeInput = document.getElementById('degree');
        const cpInput = document.getElementById('cp-count');
        const resetBtn = document.getElementById('reset-btn');
        const viewBasisBtn = document.getElementById('view-basis');
        const viewCurveBtn = document.getElementById('view-curve');

        degreeInput.addEventListener('input', (e) => {
            this.p = parseInt(e.target.value);
            document.getElementById('degree-val').textContent = this.p;
            this.rebuildKnotVector();
        });

        cpInput.addEventListener('input', (e) => {
            this.n = parseInt(e.target.value);
            document.getElementById('cp-val').textContent = this.n;
            this.points = this.initPoints();
            this.weights = new Array(this.n).fill(1.0);
            this.rebuildKnotVector();
            this.updateWeightsUI();
        });

        viewBasisBtn.addEventListener('click', () => {
            this.mode = 'basis';
            viewBasisBtn.classList.add('active');
            viewCurveBtn.classList.remove('active');
            this.render();
        });

        viewCurveBtn.addEventListener('click', () => {
            this.mode = 'curve';
            viewCurveBtn.classList.add('active');
            viewBasisBtn.classList.remove('active');
            this.render();
        });

        // Mouse Events for Interactivity
        const canvas = document.getElementById('sim-canvas');
        canvas.addEventListener('mousedown', (e) => {
            const pos = this.plot.getMousePos(e);
            // Find closest point
            let minDist = 0.05;
            this.draggingPoint = -1;
            this.points.forEach((pt, i) => {
                const d = Math.sqrt((pt.x - pos.x)**2 + (pt.y - pos.y)**2);
                if (d < minDist) {
                    minDist = d;
                    this.draggingPoint = i;
                }
            });
            this.render();
        });

        window.addEventListener('mousemove', (e) => {
            if (this.draggingPoint !== -1) {
                const pos = this.plot.getMousePos(e);
                // Clamp to [0, 1]
                this.points[this.draggingPoint].x = Math.max(0, Math.min(1, pos.x));
                this.points[this.draggingPoint].y = Math.max(0, Math.min(1, pos.y));
                this.render();
            }
        });

        window.addEventListener('mouseup', () => {
            this.draggingPoint = -1;
            this.render();
        });

        resetBtn.addEventListener('click', () => {
            this.p = 2;
            this.n = 6;
            document.getElementById('degree').value = 2;
            document.getElementById('cp-count').value = 6;
            document.getElementById('degree-val').textContent = 2;
            document.getElementById('cp-val').textContent = 6;
            this.points = this.initPoints();
            this.weights = new Array(this.n).fill(1.0);
            this.rebuildKnotVector();
            this.updateWeightsUI();
        });

        // Loop for FPS counter
        let lastTime = performance.now();
        const tick = () => {
            const now = performance.now();
            const delta = now - lastTime;
            lastTime = now;
            document.getElementById('fps').textContent = Math.round(1000 / delta);
            requestAnimationFrame(tick);
        };
        requestAnimationFrame(tick);
    }

    rebuildKnotVector() {
        this.knotVector.update(this.n, this.p);
        this.updateKnotUI();
        this.render();
    }

    updateKnotUI() {
        const editor = document.getElementById('knot-editor');
        editor.innerHTML = '';
        this.knotVector.values.forEach((val) => {
            const input = document.createElement('input');
            input.type = 'text';
            input.value = val.toFixed(2);
            input.readOnly = true;
            editor.appendChild(input);
        });
    }

    updateWeightsUI() {
        const editor = document.getElementById('weights-editor');
        editor.innerHTML = '';
        this.weights.forEach((val, i) => {
            const input = document.createElement('input');
            input.type = 'number';
            input.value = val.toFixed(1);
            input.step = 0.1;
            input.min = 0.1;
            input.addEventListener('input', (e) => {
                this.weights[i] = parseFloat(e.target.value) || 1.0;
                this.render();
            });
            editor.appendChild(input);
        });
    }

    render() {
        this.plot.clear();
        this.plot.drawGrid();
        
        if (this.mode === 'basis') {
            this.plot.drawBasisFunctions(this.n, this.p, this.knotVector.values, BasisFunctions);
        } else {
            this.plot.drawControlPolygon(this.points, this.draggingPoint);
            this.plot.drawCurve(this.p, this.knotVector.values, this.points, this.weights, Curve);
        }
    }
}

// Initialize the app
window.addEventListener('DOMContentLoaded', () => {
    new App();
});
