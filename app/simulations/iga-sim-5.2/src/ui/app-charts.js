/**
 * Phase 5.1a - Charts Manager
 * Manages sweep mapping and modal energy charts.
 */

class AppCharts {
    constructor() {
        this.sweepChart = null;
        this.energyChart = null;
        this.initSweepMap();
        this.initEnergyChart();
    }

    initSweepMap() {
        const ctx = document.getElementById('sweep-map').getContext('2d');
        this.sweepChart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Sampled Shapes (R vs L1)',
                    data: [],
                    backgroundColor: '#8b5cf6',
                    pointRadius: 5,
                    pointHoverRadius: 8
                }, {
                    label: 'Current Configuration',
                    data: [],
                    backgroundColor: '#f59e0b',
                    pointRadius: 8,
                    pointStyle: 'rectRot'
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: { 
                        title: { display: true, text: 'Radius (R)', color: '#94a3b8' },
                        grid: { color: 'rgba(255,255,255,0.05)' },
                        ticks: { color: '#64748b' }
                    },
                    y: { 
                        title: { display: true, text: 'Width (L1)', color: '#94a3b8' },
                        grid: { color: 'rgba(255,255,255,0.05)' },
                        ticks: { color: '#64748b' }
                    }
                },
                plugins: {
                    legend: { display: false }
                }
            }
        });
    }

    initEnergyChart() {
        const ctx = document.getElementById('chart-energy').getContext('2d');
        this.energyChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: [],
                datasets: [{
                    label: 'Singular Value (Log)',
                    data: [],
                    backgroundColor: 'rgba(139, 92, 246, 0.6)',
                    borderColor: '#8b5cf6',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: { 
                        type: 'logarithmic',
                        grid: { color: 'rgba(255,255,255,0.05)' },
                        ticks: { color: '#64748b' }
                    },
                    x: { 
                        grid: { display: false },
                        ticks: { color: '#64748b' }
                    }
                }
            }
        });
    }

    addSample(r, l1) {
        this.sweepChart.data.datasets[0].data.push({ x: r, y: l1 });
        this.sweepChart.update('none');
    }

    updateCurrent(r, l1) {
        this.sweepChart.data.datasets[1].data = [{ x: r, y: l1 }];
        this.sweepChart.update('none');
    }

    updateEnergy(singularValues) {
        this.energyChart.data.labels = singularValues.map((_, i) => `Mode ${i+1}`);
        this.energyChart.data.datasets[0].data = singularValues;
        this.energyChart.update();
    }

    clearSamples() {
        this.sweepChart.data.datasets[0].data = [];
        this.sweepChart.update();
    }
}

window.AppCharts = AppCharts;
