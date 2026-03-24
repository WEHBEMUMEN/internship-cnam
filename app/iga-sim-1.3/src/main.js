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

        const insertKnotSlider = document.getElementById('insert-knot');
        const insertBtn = document.getElementById('insert-btn');
        insertKnotSlider.addEventListener('input', (e) => {
            document.getElementById('knot-val').textContent = parseFloat(e.target.value).toFixed(2);
        });
        insertBtn.addEventListener('click', () => {
            const uBar = parseFloat(insertKnotSlider.value);
            this.insertKnot(uBar);
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

        document.getElementById('download-json').addEventListener('click', () => {
            const definition = {
                degree: this.p,
                controlPointsCount: this.n,
                points: this.points,
                weights: this.weights,
                knots: this.knotVector.values
            };
            const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(definition, null, 2));
            const dlAnchor = document.createElement('a');
            dlAnchor.setAttribute("href", dataStr);
            dlAnchor.setAttribute("download", "nurbs_definition.json");
            document.body.appendChild(dlAnchor);
            dlAnchor.click();
            dlAnchor.remove();
        });

        const benchmarkBtn = document.getElementById('benchmark-btn');
        let isBenchmarking = false;

        benchmarkBtn.addEventListener('click', async () => {
            if (isBenchmarking) return;
            isBenchmarking = true;
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Running...';
            
            const originalN = this.n;
            const originalP = this.p;

            // Boehm's Injection limits
            const injections = [10, 50, 100, 250, 500];
            const k_results = [];
            
            for (let K of injections) {
                this.n = 6; this.p = 2; // Baseline
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                await new Promise(r => setTimeout(r, 50));
                
                const t0 = performance.now();
                // Direct math insertion bypassing DOM for real solver metrics
                for(let i = 0; i < K; i++) {
                    const uBar = 0.1 + (Math.random() * 0.8);
                    const p = this.p;
                    const n = this.n;
                    const U = this.knotVector.values;
                    const Pw = this.points.map((pt, idx) => { const w = this.weights[idx]; return { x: pt.x * w, y: pt.y * w, w: w }; });
                    
                    let k = -1;
                    for (let j = 0; j < U.length - 1; j++) { if (uBar >= U[j] && uBar < U[j+1]) { k = j; break; } }
                    if (k === -1) continue;
                    
                    const newN = n + 1;
                    const newU = [...U]; newU.splice(k + 1, 0, uBar);
                    const newPw = [];
                    for (let j = 0; j < newN; j++) {
                        if (j <= k - p) newPw.push(Pw[j]);
                        else if (j >= k + 1) newPw.push(Pw[j - 1]);
                        else {
                            const alpha = (uBar - U[j]) / (U[j+p] - U[j]);
                            const pt0 = Pw[j-1], pt1 = Pw[j];
                            newPw.push({ x: (1 - alpha) * pt0.x + alpha * pt1.x, y: (1 - alpha) * pt0.y + alpha * pt1.y, w: (1 - alpha) * pt0.w + alpha * pt1.w });
                        }
                    }
                    this.n = newN; this.knotVector.values = newU;
                    this.points = newPw.map(pt => ({ x: pt.x / pt.w, y: pt.y / pt.w }));
                    this.weights = newPw.map(pt => pt.w);
                }
                const t1 = performance.now();
                k_results.push((t1 - t0).toFixed(1));
            }

            // Restore
            this.n = originalN;
            this.p = originalP;
            this.points = this.initPoints();
            this.weights = new Array(this.n).fill(1.0);
            this.rebuildKnotVector();
            this.updateWeightsUI();
            this.render();
            
            isBenchmarking = false;
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-gauge-high"></i> Test Knot Injection';
            
            const reportDiv = document.createElement('div');
            reportDiv.style.cssText = "position:absolute; top:50%; left:50%; transform:translate(-50%, -50%); background:var(--sidebar-bg); padding:30px; border-radius:16px; border:1px solid var(--accent); z-index:9999; box-shadow:0 20px 40px rgba(0,0,0,0.8); min-width: 400px;";
            reportDiv.innerHTML = `
                <h3 style="margin-top:0; margin-bottom: 24px; color: var(--accent); border-bottom: 1px solid var(--border); padding-bottom: 12px; font-size:1.4rem;">Boehm's Knot Insertion Matrix</h3>
                
                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test: Flooding Injection Array</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 24px;">
                    ${injections.map((K, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace;">
                            <span style="color:var(--text-secondary);">Injected knots (${K})</span>
                            <span style="color:var(--success); font-weight:bold;">${k_results[i]} ms</span>
                        </li>
                    `).join('')}
                </ul>
                <button class="btn btn-primary" style="width:100%; font-weight:bold;" onclick="this.parentElement.remove()">Close Constraints Report</button>
            `;
            document.body.appendChild(reportDiv);
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

    insertKnot(uBar) {
        const p = this.p;
        const n = this.n;
        const U = this.knotVector.values;
        const Pw = this.points.map((pt, i) => {
            const w = this.weights[i];
            return { x: pt.x * w, y: pt.y * w, w: w };
        });

        let k = -1;
        for (let i = 0; i < U.length - 1; i++) {
            if (uBar >= U[i] && uBar < U[i+1]) {
                k = i; break;
            }
        }
        if (k === -1) return; 

        const newN = n + 1;
        const newU = [...U];
        newU.splice(k + 1, 0, uBar);

        const newPw = [];
        for (let i = 0; i < newN; i++) {
            if (i <= k - p) {
                newPw.push(Pw[i]);
            } else if (i >= k + 1) {
                newPw.push(Pw[i - 1]);
            } else {
                const alpha = (uBar - U[i]) / (U[i+p] - U[i]);
                const pt0 = Pw[i-1];
                const pt1 = Pw[i];
                newPw.push({
                    x: (1 - alpha) * pt0.x + alpha * pt1.x,
                    y: (1 - alpha) * pt0.y + alpha * pt1.y,
                    w: (1 - alpha) * pt0.w + alpha * pt1.w
                });
            }
        }

        this.n = newN;
        this.knotVector.values = newU;
        this.points = newPw.map(pt => ({ x: pt.x / pt.w, y: pt.y / pt.w }));
        this.weights = newPw.map(pt => pt.w);

        document.getElementById('cp-count').max = Math.max(15, this.n + 5);
        document.getElementById('cp-count').value = this.n;
        document.getElementById('cp-val').textContent = this.n;
        this.updateKnotUI();
        this.updateWeightsUI();
        this.render();
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
