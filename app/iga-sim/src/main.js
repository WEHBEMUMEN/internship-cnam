import { KnotVector, BasisFunctions } from './engine/iga.js';
import { BasisPlot } from './ui/plot.js';

class App {
    constructor() {
        this.p = 2; // Degree
        this.n = 6; // Number of control points
        this.knotVector = new KnotVector(this.n, this.p);
        this.plot = new BasisPlot('sim-canvas');
        
        this.setupEventListeners();
        this.updateKnotUI();
        this.render();

        window.addEventListener('resize', () => this.render());
    }

    setupEventListeners() {
        const degreeInput = document.getElementById('degree');
        const cpInput = document.getElementById('cp-count');
        const resetBtn = document.getElementById('reset-btn');

        degreeInput.addEventListener('input', (e) => {
            this.p = parseInt(e.target.value);
            document.getElementById('degree-val').textContent = this.p;
            this.rebuildKnotVector();
        });

        cpInput.addEventListener('input', (e) => {
            this.n = parseInt(e.target.value);
            document.getElementById('cp-val').textContent = this.n;
            this.rebuildKnotVector();
        });

        resetBtn.addEventListener('click', () => {
            this.p = 2;
            this.n = 6;
            document.getElementById('degree').value = 2;
            document.getElementById('cp-count').value = 6;
            document.getElementById('degree-val').textContent = 2;
            document.getElementById('cp-val').textContent = 6;
            this.rebuildKnotVector();
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
        
        this.knotVector.values.forEach((val, i) => {
            const input = document.createElement('input');
            input.type = 'text';
            input.value = val.toFixed(2);
            input.readOnly = true; // For now, uniform only
            editor.appendChild(input);
        });
    }

    render() {
        this.plot.clear();
        this.plot.drawGrid();
        this.plot.drawBasisFunctions(this.n, this.p, this.knotVector.values, BasisFunctions);
    }
}

// Initialize the app
window.addEventListener('DOMContentLoaded', () => {
    new App();
});
