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
            const newP = parseInt(e.target.value);
            if (newP !== this.p) {
                this.elevateDegree(newP);
                document.getElementById('degree-val').textContent = this.p;
            }
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
            dlAnchor.setAttribute("download", "nurbs_interpolated_definition.json");
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

            // Test 1: Vary M (Fix p)
            const m_refinements = [50, 100, 250, 500, 1000];
            const fixedPForM = Math.min(this.p + 2, 6);
            const m_results = [];
            
            for (let m of m_refinements) {
                this.n = originalN; this.p = originalP;
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                await new Promise(r => setTimeout(r, 50)); 
                
                const t0 = performance.now();
                this.elevateDegree(fixedPForM, m);
                const t1 = performance.now();
                m_results.push((t1 - t0).toFixed(1));
            }

            // Test 2: Vary P (Fix M)
            const p_refinements = [2, 3, 4, 5, 6];
            const fixedM = 200;
            const p_results = [];

            for (let p of p_refinements) {
                this.n = 6; this.p = 2; 
                this.points = this.initPoints();
                this.weights = new Array(this.n).fill(1.0);
                this.knotVector.update(this.n, this.p);
                
                await new Promise(r => setTimeout(r, 50));
                
                const t0 = performance.now();
                this.elevateDegree(p, fixedM);
                const t1 = performance.now();
                p_results.push((t1 - t0).toFixed(1));
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
            benchmarkBtn.innerHTML = '<i class="fa-solid fa-gauge-high"></i> Test Matrix Solver';
            
            const reportDiv = document.createElement('div');
            reportDiv.style.cssText = "position:absolute; top:50%; left:50%; transform:translate(-50%, -50%); background:var(--sidebar-bg); padding:30px; border-radius:16px; border:1px solid var(--accent); z-index:9999; box-shadow:0 20px 40px rgba(0,0,0,0.8); min-width: 400px;";
            reportDiv.innerHTML = `
                <h3 style="margin-top:0; margin-bottom: 24px; color: var(--accent); border-bottom: 1px solid var(--border); padding-bottom: 12px; font-size:1.4rem;">Linear Algebra Benchmark</h3>
                
                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 1: Scaling Samples (Fixed p=${fixedPForM})</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 24px;">
                    ${m_refinements.map((m, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace;">
                            <span style="color:var(--text-secondary);">Samples (M=${m})</span>
                            <span style="color:var(--success); font-weight:bold;">${m_results[i]} ms</span>
                        </li>
                    `).join('')}
                </ul>

                <h4 style="margin:0 0 12px 0; color:var(--text-secondary); font-size:0.85rem; text-transform:uppercase; letter-spacing:1px;">Test 2: Scaling Degree (Fixed M=${fixedM})</h4>
                <ul style="list-style:none; padding:0; margin-bottom: 30px;">
                    ${p_refinements.map((p, i) => `
                        <li style="display:flex; justify-content:space-between; margin-bottom:8px; font-family:'JetBrains Mono',monospace;">
                            <span style="color:var(--text-secondary);">Degree (p=${p})</span>
                            <span style="color:var(--success); font-weight:bold;">${p_results[i]} ms</span>
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

    elevateDegree(newP, M = 150) {
        // Collect samples from current curve
        const samples = [];
        for(let i=0; i<=M; i++) {
            const xi = i/M;
            samples.push(Curve.evaluate(xi, this.p, this.knotVector.values, this.points, this.weights));
        }

        // Increase control points to maintain shape accuracy
        const newN = this.n + Math.abs(newP - this.p) * 2; 
        const newKnotVector = new KnotVector(newN, newP);
        const newU = newKnotVector.values;
        const newWeights = new Array(newN).fill(1.0);

        // Build A matrix (M+1 x newN)
        const A = [];
        for(let i=0; i<=M; i++) {
            const xi = i/M;
            A.push(NURBS.evaluateAll(newN, newP, xi, newU, newWeights));
        }

        // Build ATA (newN x newN) and ATB (newN x 2)
        const ATA = Array(newN).fill(0).map(()=>Array(newN).fill(0));
        const ATBx = Array(newN).fill(0);
        const ATBy = Array(newN).fill(0);

        for(let i=0; i<newN; i++) {
            for(let j=0; j<newN; j++) {
                let sum = 0;
                for(let k=0; k<=M; k++) { sum += A[k][i] * A[k][j]; }
                ATA[i][j] = sum;
            }
            let sumX = 0, sumY = 0;
            for(let k=0; k<=M; k++) {
                sumX += A[k][i] * samples[k].x;
                sumY += A[k][i] * samples[k].y;
            }
            ATBx[i] = sumX;
            ATBy[i] = sumY;
        }

        // Solve ATA * X = ATBx
        const solve = (matA, matB) => {
            let n = matA.length;
            let A = matA.map(r => [...r]);
            let B = [...matB];
            for (let i = 0; i < n; i++) {
                let maxEl = Math.abs(A[i][i]), maxRow = i;
                for (let k = i + 1; k < n; k++) {
                    if (Math.abs(A[k][i]) > maxEl) { maxEl = Math.abs(A[k][i]); maxRow = k; }
                }
                for (let k = i; k < n; k++) {
                    let tmp = A[maxRow][k]; A[maxRow][k] = A[i][k]; A[i][k] = tmp;
                }
                let tmpB = B[maxRow]; B[maxRow] = B[i]; B[i] = tmpB;
                
                for (let k = i + 1; k < n; k++) {
                    if (A[i][i] === 0) continue;
                    let c = -A[k][i] / A[i][i];
                    for (let j = i; j < n; j++) {
                        if (i === j) A[k][j] = 0; else A[k][j] += c * A[i][j];
                    }
                    B[k] += c * B[i];
                }
            }
            let x = new Array(n).fill(0);
            for (let i = n - 1; i > -1; i--) {
                if (A[i][i] === 0) { x[i] = 0; continue; }
                let sum = 0;
                for (let k = i + 1; k < n; k++) { sum += A[i][k] * x[k]; }
                x[i] = (B[i] - sum) / A[i][i];
            }
            return x;
        };

        const newX = solve(ATA, ATBx);
        const newY = solve(ATA, ATBy);

        const newPoints = [];
        for(let i=0; i<newN; i++) { newPoints.push({x: newX[i], y: newY[i]}); }
        if (newN > 0 && samples.length > 0) {
           newPoints[0] = {x: samples[0].x, y: samples[0].y};
           newPoints[newN - 1] = {x: samples[M].x, y: samples[M].y};
        }

        this.points = newPoints;
        this.n = newN;
        this.weights = newWeights;
        this.p = newP;
        this.knotVector = newKnotVector;
        
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
