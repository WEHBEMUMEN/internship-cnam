/**
 * Phase 3.4b: U-DEIM Benchmark — UI & Charts
 */

UDEIMBenchmarkApp.prototype.initUI = function() {
    document.getElementById('btn-train').onclick = () => this.trainAll();
    document.getElementById('btn-compare').onclick = () => this.runComparison();
    document.getElementById('btn-audit').onclick = () => {
        if (this.isTrained) this.deimEngine.auditMath(this.solverFOM, this.romEngine, this.patch, this.snapDisp);
    };
    
    const btnExp = document.getElementById('btn-explorer');
    btnExp.onclick = () => {
        if (!this.isTrained) return;
        document.getElementById('tab-explorer').style.display = 'block';
        document.getElementById('tab-explorer').click();
    };

    document.getElementById('btn-exp-prev').onclick = () => {
        if (this.explorerStep > 1) { this.explorerStep--; this.runDEIMExplorer(this.explorerStep); }
    };
    document.getElementById('btn-exp-next').onclick = () => {
        if (this.explorerStep < this.deimEngine.history.length) { this.explorerStep++; this.runDEIMExplorer(this.explorerStep); }
    };

    document.querySelectorAll('[data-method]').forEach(btn => {
        btn.onclick = () => {
            document.querySelectorAll('[data-method]').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            this.method = btn.dataset.method;
            if (this.method === 'fom' || this.isTrained) this.updatePhysics();
        };
    });

    document.getElementById('input-load').oninput = e => {
        this.loadMag = parseInt(e.target.value);
        document.getElementById('load-val').textContent = this.loadMag;
        this._scheduleUpdate();
    };

    document.getElementById('input-k').oninput = e => {
        this.k = parseInt(e.target.value);
        document.getElementById('k-val').textContent = this.k;
        if (this.isTrained) {
            this.romEngine.computePOD(this.k, this.patch, this.getBCs());

            if (this.snapDisp) this.deimEngine.train(this.solverFOM, this.romEngine, this.patch, this.snapDisp, this.deimM, this.k).then(() => {
                this._scheduleUpdate();
                this.deimEngine.auditMath(this.solverFOM, this.romEngine, this.patch, this.snapDisp);
            });

        }
    };

    document.getElementById('input-m').oninput = e => {
        this.deimM = parseInt(e.target.value);
        document.getElementById('m-val').textContent = this.deimM;
        if (this.isTrained && this.snapDisp) {
            this.deimEngine.train(this.solverFOM, this.romEngine, this.patch, this.snapDisp, this.deimM, this.k).then(() => {
                this._scheduleUpdate();
                this.deimEngine.auditMath(this.solverFOM, this.romEngine, this.patch, this.snapDisp);
            });
        }

    };

    document.getElementById('input-mesh').onchange = e => {
        this.meshLevel = parseInt(e.target.value);
        document.getElementById('mesh-val').textContent = this.meshLevel;
        this.loadBenchmark();
    };

    document.querySelectorAll('[data-load-type]').forEach(btn => {
        btn.onclick = () => {
            document.querySelectorAll('[data-load-type]').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            this.loadType = btn.dataset.loadType;
            this.updatePhysics();
        };
    });
};

UDEIMBenchmarkApp.prototype._scheduleUpdate = function() {
    if (this._timer) clearTimeout(this._timer);
    this._timer = setTimeout(() => {
        this.updatePhysics();
        if ((this.currentChartTab === 'fd' || this.currentChartTab === 'error') && this.isTrained) {
            this.runFDCurves();
        }
    }, 100);
};


UDEIMBenchmarkApp.prototype._updateStats = function(method, meta) {
    document.getElementById('method-name').textContent = METHODS[method].label;
    document.getElementById('time-val').textContent = meta.time.toFixed(1);
    const speedup = this.lastFomTime ? (this.lastFomTime / meta.time).toFixed(1) : '—';
    document.getElementById('speedup-val').textContent = method === 'fom' ? '1.0' : speedup;
    document.getElementById('sampled-val').textContent = meta.sampled || 'All';
    document.getElementById('tip-val').textContent = meta.tipDisp ? Math.abs(meta.tipDisp).toFixed(3) : '—';
    document.getElementById('error-val').textContent = meta.error !== undefined ? (meta.error*100).toFixed(2)+'%' : '—';
};

UDEIMBenchmarkApp.prototype.initCharts = function() {
    this.sparkline = new Sparkline(document.getElementById('sparkline-canvas'));
    this.speedupChart = new Chart(document.getElementById('chart-speedup'), {
        type: 'bar',
        data: { labels: [], datasets: [{ label: 'Speedup (×)', data: [], backgroundColor: [] }] },
        options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{display:false}},
            scales:{y:{beginAtZero:true, title:{display:true, text:'Speedup (×)', font:{size:10}}}, x:{ticks:{font:{size:9}}}}}
    });
    this.fdChart = new Chart(document.getElementById('chart-fd'), {
        type: 'line', data: { labels: [], datasets: [] },
        options: { responsive:true, maintainAspectRatio:false,
            plugins:{legend:{position:'top', labels:{font:{size:9}}}},
            scales:{x:{title:{display:true, text:'Load F', font:{size:10}}}, y:{title:{display:true, text:'Tip Displacement', font:{size:10}}}}}
    });
    document.querySelectorAll('.chart-tab').forEach(t => {
        t.onclick = () => {
            document.querySelectorAll('.chart-tab').forEach(b => b.classList.remove('active'));
            t.classList.add('active');
            document.getElementById('chart-speedup-wrap').classList.toggle('hidden', t.dataset.chart !== 'speedup');
            document.getElementById('chart-fd-wrap').classList.toggle('hidden', t.dataset.chart !== 'fd');
            document.getElementById('chart-error-wrap').classList.toggle('hidden', t.dataset.chart !== 'error');
            document.getElementById('explorer-wrap').classList.toggle('hidden', t.dataset.chart !== 'explorer');
            
            this.currentChartTab = t.dataset.chart;
            if ((t.dataset.chart === 'fd' || t.dataset.chart === 'error') && this.isTrained) this.runFDCurves();

            if (t.dataset.chart === 'explorer') {
                this.isExplorerActive = true;
                this.runDEIMExplorer(this.explorerStep);
            } else {
                this.isExplorerActive = false;
                this.sensorSpheres.clear();
                if (this.activeElementMesh) this.activeElementMesh.visible = false;
                if (this.lastResult) this.updateMesh(this.lastResult.u);
                this._render();
            }
        };
    });

    this.errorChart = new Chart(document.getElementById('chart-error'), {
        type: 'line', data: { labels: [], datasets: [] },
        options: { responsive:true, maintainAspectRatio:false,
            plugins:{legend:{position:'top', labels:{font:{size:9}}}},
            scales:{x:{title:{display:true, text:'Load F', font:{size:10}}}, y:{type: 'logarithmic', title:{display:true, text:'L2 Rel Error (%)', font:{size:10}}, ticks:{callback: (val) => val.toFixed(4) + '%' }}}}


    });

    this.residualChart = new Chart(document.getElementById('chart-residual'), {
        type: 'line',
        data: { labels: [], datasets: [
            { label: 'Residual', data: [], borderColor: '#f59e0b', backgroundColor: 'rgba(245,158,11,0.2)', borderWidth: 2, tension: 0.1, fill: true, pointRadius: 0 },
            { label: 'Selected Node', data: [], pointBackgroundColor: '#ef4444', pointBorderColor: '#fff', pointRadius: 5, pointHoverRadius: 7, showLine: false }
        ]},
        options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{display:false}, tooltip:{enabled:false}}, animation: false,
            scales:{x:{display:false}, y:{type:'logarithmic', title:{display:true, text:'Magnitude', font:{size:9}}, ticks:{font:{size:8}}}}
        }
    });
};

UDEIMBenchmarkApp.prototype.updateSpeedupChart = function(data) {
    const labels = [], values = [], bgColors = [];
    const fomTime = data.fom?.time || 1;
    for (const [m, meta] of Object.entries(data)) {
        labels.push(METHODS[m].label);
        values.push(parseFloat((fomTime / meta.time).toFixed(1)));
        bgColors.push(METHODS[m].color + '99');
    }
    this.speedupChart.data.labels = labels;
    this.speedupChart.data.datasets[0].data = values;
    this.speedupChart.data.datasets[0].backgroundColor = bgColors;
    this.speedupChart.update();
};
