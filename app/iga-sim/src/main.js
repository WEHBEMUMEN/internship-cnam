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
        
        // Performance & Resolution State
        this.resolution = 500;
        this.showPoints = false;
        this.evalTime = 0;
        
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

        const resSlider = document.getElementById('input-res');
        if (resSlider) {
            resSlider.addEventListener('input', (e) => {
                this.resolution = parseInt(e.target.value);
                document.getElementById('res-val').textContent = this.resolution;
                this.render();
            });
        }

        const showPointsCheck = document.getElementById('check-show-points');
        if (showPointsCheck) {
            showPointsCheck.addEventListener('change', (e) => {
                this.showPoints = e.target.checked;
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

        document.getElementById('import-json-btn').addEventListener('click', () => document.getElementById('import-json-file').click());
        document.getElementById('import-json-file').addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (!file) return;
            const reader = new FileReader();
            reader.onload = (event) => {
                try {
                    const data = JSON.parse(event.target.result);
                    this.p = data.degree || data.p;
                    this.n = data.controlPointsCount || data.n;
                    this.points = data.points;
                    this.weights = data.weights;
                    this.knotVector.update(this.n, this.p);
                    this.knotVector.values = data.knots;
                    document.getElementById('degree').value = this.p;
                    document.getElementById('degree-val').textContent = this.p;
                    document.getElementById('cp-count').value = this.n;
                    document.getElementById('cp-val').textContent = this.n;
                    this.updateKnotUI();
                    this.updateWeightsUI();
                    this.render();
                } catch (err) { alert("Invalid Matrix Data Frame format."); }
            };
            reader.readAsText(file);
        });

        const benchmarkBtn = document.getElementById('benchmark-btn');
        let isBenchmarking = false;

        benchmarkBtn.addEventListener('click', async () => {
            if (isBenchmarking) return;
            isBenchmarking = true;
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Running...';
            
            const originalN = this.n;
            const originalP = this.p;
            const originalRes = this.resolution;
            
            // --- Test 1: Scaling Resolution (Evaluation Throughput) ---
            const res_refinements = [10, 100, 500, 1000, 5000];
            const res_results = [];
            
            for (let res of res_refinements) {
                this.resolution = res;
                let frames = 0;
                let totalEval = 0;
                
                const startTime = performance.now();
                await new Promise(resolve => {
                    const loop = () => {
                        if (frames >= 30) { resolve(); return; }
                        const t0 = performance.now();
                        this.render(); 
                        totalEval += (performance.now() - t0);
                        frames++;
                        requestAnimationFrame(loop);
                    };
                    requestAnimationFrame(loop);
                });
                
                const endTime = performance.now();
                res_results.push({
                    fps: Math.floor(30 / ((endTime - startTime) / 1000)),
                    avgCpu: (totalEval / 30).toFixed(3)
                });
            }
            
            // --- Test 2: Scaling Control Points at High Res ---
            const n_refinements = [6, 20, 50, 100];
            const n_results = [];
            this.resolution = 1000;

            for (let n of n_refinements) {
                this.n = n;
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                let frames = 0;
                let totalEval = 0;
                const startTime = performance.now();
                await new Promise(resolve => {
                    const loop = () => {
                        if (frames >= 30) { resolve(); return; }
                        const t0 = performance.now();
                        this.render(); 
                        totalEval += (performance.now() - t0);
                        frames++;
                        requestAnimationFrame(loop);
                    };
                    requestAnimationFrame(loop);
                });
                const endTime = performance.now();
                n_results.push({
                    fps: Math.floor(30 / ((endTime - startTime) / 1000)),
                    avgCpu: (totalEval / 30).toFixed(3)
                });
            }

            // Restore Safe State
            this.n = originalN;
            this.p = originalP;
            this.resolution = originalRes;
            this.points = this.initPoints();
            this.weights = new Array(this.n).fill(1.0);
            this.rebuildKnotVector();
            this.updateWeightsUI();
            this.render();
            
            isBenchmarking = false;
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-gauge-high"></i> Benchmark Test';
            
            // Generate Unified Multi-Variable Report
            const reportDiv = document.createElement('div');
            reportDiv.style.cssText = "position:absolute; top:50%; left:50%; transform:translate(-50%, -50%); background:var(--sidebar-bg); padding:30px; border-radius:16px; border:1px solid var(--accent); z-index:9999; box-shadow:0 20px 40px rgba(0,0,0,0.8); min-width: 450px;";
            reportDiv.innerHTML = `
                <h3 style="margin-top:0; margin-bottom: 24px; color: var(--accent); border-bottom: 1px solid var(--border); padding-bottom: 12px; font-size:1.4rem;">Performance Architecture Report</h3>
                
                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 1: Sampling Scaling (Fixed n, p)</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 24px;">
                    ${res_refinements.map((res, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace; font-size:0.85rem;">
                            <span style="color:var(--text-secondary);">Pts: ${res}</span>
                            <span>
                                <span style="color:var(--accent); margin-right:15px;">${res_results[i].avgCpu}ms</span>
                                <span style="color:var(--success); font-weight:bold;">${res_results[i].fps} FPS</span>
                            </span>
                        </li>
                    `).join('')}
                </ul>

                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 2: Complexity Scaling (at 1000 Pts)</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 30px;">
                    ${n_refinements.map((n, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace; font-size:0.85rem;">
                            <span style="color:var(--text-secondary);">Nodes: ${n}</span>
                            <span>
                                <span style="color:var(--accent); margin-right:15px;">${n_results[i].avgCpu}ms</span>
                                <span style="color:var(--success); font-weight:bold;">${n_results[i].fps} FPS</span>
                            </span>
                        </li>
                    `).join('')}
                </ul>

                <button class="btn btn-primary" style="width:100%; font-weight:bold;" onclick="this.parentElement.remove()">Close Performance Report</button>
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
        const t0 = performance.now();
        this.plot.clear();
        this.plot.drawGrid();
        
        if (this.mode === 'basis') {
            this.plot.drawBasisFunctions(this.n, this.p, this.knotVector.values, BasisFunctions, this.resolution);
        } else {
            this.plot.drawControlPolygon(this.points, this.draggingPoint);
            this.plot.drawCurve(this.p, this.knotVector.values, this.points, this.weights, Curve, this.resolution, this.showPoints);
        }
        const t1 = performance.now();
        this.evalTime = t1 - t0;
        
        const cpuEl = document.getElementById('cpu-time');
        if (cpuEl) cpuEl.textContent = this.evalTime.toFixed(3);
        const ptsEl = document.getElementById('pts-count');
        if (ptsEl) ptsEl.textContent = this.resolution;
    }
}

// Initialize the app
window.addEventListener('DOMContentLoaded', () => {
    new App();
});
