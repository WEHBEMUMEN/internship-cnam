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
}
