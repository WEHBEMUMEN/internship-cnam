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
        
        // Horizontal lines
        for (let i = 0; i <= 4; i++) {
            const y = height * (1 - i/4) * 0.8 + height * 0.1;
            this.ctx.beginPath();
            this.ctx.moveTo(width * 0.1, y);
            this.ctx.lineTo(width * 0.9, y);
            this.ctx.stroke();
        }
        
        // Vertical axis
        this.ctx.beginPath();
        this.ctx.moveTo(width * 0.1, height * 0.1);
        this.ctx.lineTo(width * 0.1, height * 0.9);
        this.ctx.stroke();
        
        // Horizontal axis
        this.ctx.beginPath();
        this.ctx.moveTo(width * 0.1, height * 0.9);
        this.ctx.lineTo(width * 0.9, height * 0.9);
        this.ctx.stroke();
    }

    /**
     * Draws all basis functions
     */
    drawBasisFunctions(n, p, knots, BasisEngine) {
        const { width, height } = this.canvas;
        const colors = [
            '#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', 
            '#ec4899', '#06b6d4', '#f97316', '#22c55e', '#a855f7'
        ];

        for (let i = 0; i < n; i++) {
            const color = colors[i % colors.length];
            this.ctx.strokeStyle = color;
            this.ctx.lineWidth = 3;
            this.ctx.lineJoin = 'round';
            this.ctx.beginPath();

            let first = true;
            // Sampling across Xi interval [0, 1]
            for (let x = 0; x <= 1; x += 0.005) {
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
            this.ctx.lineTo(width * 0.9, height * 0.9); // Close to bottom right
            this.ctx.lineTo(width * 0.1, height * 0.9); // Back to bottom left
            this.ctx.fill();
        }
    }

    /**
     * Draws the physical NURBS curve C(xi)
     */
    drawCurve(p, knots, points, weights, CurveEngine) {
        const { width, height } = this.canvas;
        this.ctx.strokeStyle = '#ffffff';
        this.ctx.lineWidth = 4;
        this.ctx.lineJoin = 'round';
        this.ctx.setLineDash([]); // Ensure solid line

        this.ctx.beginPath();
        let first = true;

        for (let xi = 0; xi <= 1; xi += 0.002) {
            const pt = CurveEngine.evaluate(xi, p, knots, points, weights);
            // In 1D, pt.x is parameter, pt.y is the mapped value
            // But here the user wants "interactive points" to draw a line.
            // We will treat points[i].x and points[i].y as actual screen coords (normalized 0-1)
            const px = width * (0.1 + pt.x * 0.8);
            const py = height * (0.1 + pt.y * 0.8);

            if (first) {
                this.ctx.moveTo(px, py);
                first = false;
            } else {
                this.ctx.lineTo(px, py);
            }
        }
        this.ctx.stroke();
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
