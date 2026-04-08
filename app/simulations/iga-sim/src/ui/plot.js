/**
 * IGA Plot Renderer
 * High-performance Canvas rendering for basis functions
 */

export class BasisPlot {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        this.resize();
        window.addEventListener('resize', () => this.resize());
    }

    resize() {
        const rect = this.canvas.parentElement.getBoundingClientRect();
        this.canvas.width = rect.width;
        this.canvas.height = rect.height;
    }

    clear() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    }

    /**
     * Draws the coordinate system
     */
    drawGrid() {
        const { width, height } = this.canvas;
        this.ctx.strokeStyle = '#334155';
        this.ctx.lineWidth = 1;
        this.ctx.font = '12px "JetBrains Mono", monospace';
        this.ctx.fillStyle = '#94a3b8';
        this.ctx.textAlign = 'center';
        
        // Horizontal lines (Grid)
        for (let i = 0; i <= 4; i++) {
            const y = height * (1 - i/4) * 0.8 + height * 0.1;
            this.ctx.beginPath();
            this.ctx.moveTo(width * 0.1, y);
            this.ctx.lineTo(width * 0.9, y);
            this.ctx.stroke();
            
            // Y-Axis Ticks
            const val = (i/4).toFixed(2);
            if (i % 2 === 0) { // Labels for 0, 0.5, 1.0
                this.ctx.fillText(val, width * 0.07, y + 4);
            }
        }
        
        // Vertical axis (Y)
        this.ctx.beginPath();
        this.ctx.lineWidth = 2;
        this.ctx.strokeStyle = '#475569';
        this.ctx.moveTo(width * 0.1, height * 0.9);
        this.ctx.lineTo(width * 0.1, height * 0.05);
        this.ctx.stroke();
        
        // Y arrowhead
        this.ctx.beginPath();
        this.ctx.moveTo(width * 0.1 - 5, height * 0.07);
        this.ctx.lineTo(width * 0.1, height * 0.05);
        this.ctx.lineTo(width * 0.1 + 5, height * 0.07);
        this.ctx.stroke();

        // Y label
        const yLabel = document.getElementById('view-basis')?.classList.contains('active') ? 'N(ξ)' : 'y';
        this.ctx.fillText(yLabel, width * 0.1, height * 0.03);
        
        // Horizontal axis (X)
        this.ctx.beginPath();
        this.ctx.moveTo(width * 0.1, height * 0.9);
        this.ctx.lineTo(width * 0.95, height * 0.9);
        this.ctx.stroke();
        
        // X arrowhead
        this.ctx.beginPath();
        this.ctx.moveTo(width * 0.93, height * 0.9 - 5);
        this.ctx.lineTo(width * 0.95, height * 0.9);
        this.ctx.lineTo(width * 0.93, height * 0.9 + 5);
        this.ctx.stroke();

        // X label
        this.ctx.fillText('ξ', width * 0.97, height * 0.9 + 4);

        // X ticks (0, 0.5, 1.0)
        [0, 0.5, 1.0].forEach(val => {
            const tx = width * 0.1 + val * width * 0.8;
            this.ctx.fillText(val.toFixed(1), tx, height * 0.93);
        });
    }

    /**
     * Draws all basis functions with variable resolution
     */
    drawBasisFunctions(n, p, knots, BasisEngine, resolution = 200) {
        const { width, height } = this.canvas;
        const colors = [
            '#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', 
            '#ec4899', '#06b6d4', '#f97316', '#22c55e', '#a855f7'
        ];

        const step = 1.0 / (resolution - 1);

        for (let i = 0; i < n; i++) {
            const color = colors[i % colors.length];
            this.ctx.strokeStyle = color;
            this.ctx.lineWidth = 3;
            this.ctx.lineJoin = 'round';
            this.ctx.beginPath();

            let first = true;
            for (let j = 0; j < resolution; j++) {
                const x = j * step;
                const val = BasisEngine.evaluate(i, p, x, knots);
                const px = width * 0.1 + x * width * 0.8;
                const py = height * 0.9 - val * height * 0.7;
                
                if (first) {
                    this.ctx.moveTo(px, py);
                    first = false;
                } else {
                    this.ctx.lineTo(px, py);
                }
            }
            this.ctx.stroke();
            
            // Draw area under curve with transparency
            this.ctx.fillStyle = color + '22';
            this.ctx.lineTo(width * 0.9, height * 0.9);
            this.ctx.lineTo(width * 0.1, height * 0.9);
            this.ctx.fill();
        }
    }

    /**
     * Draws the physical NURBS curve C(xi) with variable resolution
     */
    drawCurve(p, knots, points, weights, CurveEngine, resolution = 500, showPoints = false) {
        const { width, height } = this.canvas;
        this.ctx.strokeStyle = '#ffffff';
        this.ctx.lineWidth = 4;
        this.ctx.lineJoin = 'round';
        this.ctx.setLineDash([]);

        this.ctx.beginPath();
        let first = true;

        const step = 1.0 / (resolution - 1);
        const evaluatedPoints = [];

        for (let i = 0; i < resolution; i++) {
            const xi = i * step;
            const pt = CurveEngine.evaluate(xi, p, knots, points, weights);
            const px = width * (0.1 + pt.x * 0.8);
            const py = height * (0.1 + pt.y * 0.8);
            
            evaluatedPoints.push({x: px, y: py});

            if (first) {
                this.ctx.moveTo(px, py);
                first = false;
            } else {
                this.ctx.lineTo(px, py);
            }
        }
        this.ctx.stroke();

        // Optional: Draw evaluation points as dots
        if (showPoints) {
            this.ctx.fillStyle = '#3b82f6';
            evaluatedPoints.forEach(pt => {
                this.ctx.beginPath();
                this.ctx.arc(pt.x, pt.y, 3, 0, Math.PI * 2);
                this.ctx.fill();
            });
        }
    }

    /**
     * Draws the control polygon and points
     */
    drawControlPolygon(points, activeIndex = -1) {
        const { width, height } = this.canvas;
        
        // Draw lines between points
        this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
        this.ctx.setLineDash([5, 5]);
        this.ctx.lineWidth = 1;
        this.ctx.beginPath();
        points.forEach((pt, i) => {
            const px = width * (0.1 + pt.x * 0.8);
            const py = height * (0.1 + pt.y * 0.8);
            if (i === 0) this.ctx.moveTo(px, py);
            else this.ctx.lineTo(px, py);
        });
        this.ctx.stroke();
        this.ctx.setLineDash([]);

        // Draw points
        points.forEach((pt, i) => {
            const px = width * (0.1 + pt.x * 0.8);
            const py = height * (0.1 + pt.y * 0.8);
            
            this.ctx.fillStyle = i === activeIndex ? '#3b82f6' : '#94a3b8';
            this.ctx.beginPath();
            this.ctx.arc(px, py, i === activeIndex ? 8 : 6, 0, Math.PI * 2);
            this.ctx.fill();
            
            // Halo for active point
            if (i === activeIndex) {
                this.ctx.strokeStyle = 'rgba(59, 130, 246, 0.4)';
                this.ctx.lineWidth = 4;
                this.ctx.stroke();
            }
        });
    }

    /**
     * Gets normalized coordinates from mouse event
     */
    getMousePos(e) {
        const rect = this.canvas.getBoundingClientRect();
        const x = (e.clientX - rect.left) / rect.width;
        const y = (e.clientY - rect.top) / rect.height;
        
        // Normalize back to the [0, 1] range inside our 0.1-0.9 padding
        return {
            x: (x - 0.1) / 0.8,
            y: (y - 0.1) / 0.8
        };
    }
}
