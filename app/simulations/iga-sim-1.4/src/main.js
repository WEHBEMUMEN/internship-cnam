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
        this.evalPoint = 0.5;
        
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
        const evalPointInput = document.getElementById('eval-point');

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

        if (evalPointInput) {
            evalPointInput.addEventListener('input', (e) => {
                this.evalPoint = parseFloat(e.target.value);
                document.getElementById('eval-val').textContent = this.evalPoint.toFixed(2);
                this.render();
            });
        }

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
            
            // --- Test 1: Vary Control Points (Fix Degree at Current) ---
            const n_refinements = [10, 50, 100, 200, 500];
            const n_results = [];
            
            for (let n of n_refinements) {
                this.n = n;
                this.p = originalP;
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                let frames = 0;
                let angle = 0;
                
                const startTime = performance.now();
                await new Promise(resolve => {
                    const loop = () => {
                        if (frames >= 60) { resolve(); return; }
                        this.weights = this.weights.map((w, i) => 1.0 + 0.8 * Math.sin(angle + i));
                        this.points.forEach((pt, i) => { pt.y = 0.5 + 0.3 * Math.sin(angle * 2 + i); });
                        this.render(); 
                        angle += 0.2;
                        frames++;
                        requestAnimationFrame(loop);
                    };
                    requestAnimationFrame(loop);
                });
                
                const endTime = performance.now();
                n_results.push(Math.floor(60 / ((endTime - startTime) / 1000)));
            }
            
            // --- Test 2: Vary Degree (Fix Control Points at safe constant N=50) ---
            const p_refinements = [2, 3, 4, 5, 6]; 
            const fixedN = Math.max(originalN, 15);
            const p_results = [];

            for (let p of p_refinements) {
                this.p = p;
                this.n = fixedN; 
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                let frames = 0;
                let angle = 0;
                
                const startTime = performance.now();
                await new Promise(resolve => {
                    const loop = () => {
                        if (frames >= 60) { resolve(); return; }
                        this.weights = this.weights.map((w, i) => 1.0 + 0.8 * Math.sin(angle + i));
                        this.points.forEach((pt, i) => { pt.y = 0.5 + 0.3 * Math.sin(angle * 2 + i); });
                        this.render(); 
                        angle += 0.2;
                        frames++;
                        requestAnimationFrame(loop);
                    };
                    requestAnimationFrame(loop);
                });
                
                const endTime = performance.now();
                p_results.push(Math.floor(60 / ((endTime - startTime) / 1000)));
            }

            // Restore Safe State
            this.n = originalN;
            this.p = originalP;
            this.points = this.initPoints();
            this.weights = new Array(this.n).fill(1.0);
            this.rebuildKnotVector();
            this.updateWeightsUI();
            this.render();
            
            isBenchmarking = false;
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-gauge-high"></i> Benchmark Test';
            
            // Generate Unified Multi-Variable Report
            const reportDiv = document.createElement('div');
            reportDiv.style.cssText = "position:absolute; top:50%; left:50%; transform:translate(-50%, -50%); background:var(--sidebar-bg); padding:30px; border-radius:16px; border:1px solid var(--accent); z-index:9999; box-shadow:0 20px 40px rgba(0,0,0,0.8); min-width: 360px;";
            reportDiv.innerHTML = `
                <h3 style="margin-top:0; margin-bottom: 24px; color: var(--accent); border-bottom: 1px solid var(--border); padding-bottom: 12px; font-size:1.4rem;">Isolated Architecture Benchmark</h3>
                
                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 1: Scaling Nodes (Fixed p=${originalP})</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 24px;">
                    ${n_refinements.map((n, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace;">
                            <span style="color:var(--text-secondary);">Nodes (n=${n})</span>
                            <span style="color:var(--success); font-weight:bold;">${n_results[i]} FPS</span>
                        </li>
                    `).join('')}
                </ul>

                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 2: Scaling Degree (Fixed n=${fixedN})</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 30px;">
                    ${p_refinements.map((p, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace;">
                            <span style="color:var(--text-secondary);">Degree (p=${p})</span>
                            <span style="color:var(--success); font-weight:bold;">${p_results[i]} FPS</span>
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
        if (window.perfMonitor) window.perfMonitor.startMeasure();
        this.plot.clear();
        this.plot.drawGrid();
        
        if (this.mode === 'basis') {
            this.plot.drawBasisFunctions(this.n, this.p, this.knotVector.values, BasisFunctions);
        } else {
            this.plot.drawControlPolygon(this.points, this.draggingPoint);
            this.plot.drawCurve(this.p, this.knotVector.values, this.points, this.weights, Curve);
            this.plot.drawCalculusVectors(this.evalPoint, this.p, this.knotVector.values, this.points, this.weights, Curve);
        }
        if (window.perfMonitor) window.perfMonitor.endMeasure();
    }
}

// Initialize the app
window.addEventListener('DOMContentLoaded', () => {
    new App();
});
