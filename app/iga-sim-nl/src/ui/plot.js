/**
 * IGA Plot Renderer - Simplified & Robust for NL App
 */
export class BasisPlot {
    constructor(canvasOrId) {
        this.canvas = typeof canvasOrId === 'string' ? document.getElementById(canvasOrId) : canvasOrId;
        this.ctx = this.canvas.getContext('2d');
        this.camera = { x: 0, y: 0, zoom: 1.0 };
        this.mousePos = { x: -1000, y: -1000, worldX: 0, worldY: 0 };
        this.hoverValue = null;
        
        this.setupEvents();
        this.resize();
        new ResizeObserver(() => this.resize()).observe(this.canvas.parentElement);
    }

    setupEvents() {
        this.canvas.addEventListener('mousemove', (e) => {
            const rect = this.canvas.getBoundingClientRect();
            this.mousePos.x = e.clientX - rect.left;
            this.mousePos.y = e.clientY - rect.top;
            const wPos = this.screenToWorld(this.mousePos.x, this.mousePos.y);
            this.mousePos.worldX = wPos.x;
            this.mousePos.worldY = wPos.y;
            
            const readout = document.getElementById('coord-readout');
            if (readout) readout.textContent = `ξ: ${this.mousePos.worldX.toFixed(3)} | y: ${this.mousePos.worldY.toFixed(3)}`;
        });
        
        document.getElementById('btn-zoom-in')?.addEventListener('click', () => { this.camera.zoom *= 1.2; });
        document.getElementById('btn-zoom-out')?.addEventListener('click', () => { this.camera.zoom /= 1.2; });
        document.getElementById('btn-zoom-reset')?.addEventListener('click', () => { this.camera = { x: 0, y: 0, zoom: 1.0 }; });
    }

    resize() {
        if (!this.canvas) return;
        const rect = this.canvas.parentElement.getBoundingClientRect();
        this.canvas.width = rect.width;
        this.canvas.height = rect.height;
    }

    worldToScreen(wx, wy) {
        const { width, height } = this.canvas;
        const scaleX = width * 0.8 * this.camera.zoom;
        const scaleY = height * 0.8 * this.camera.zoom;
        
        // We want (0.5, 0.5) to be screen center
        // worldX 0 to 1 => screenX 0.1w to 0.9w
        const sx = width / 2 + (wx - 0.5 + this.camera.x) * scaleX;
        const sy = height / 2 - (wy - 0.5 + this.camera.y) * scaleY;
        return { x: sx, y: sy };
    }

    screenToWorld(sx, sy) {
        const { width, height } = this.canvas;
        const scaleX = width * 0.8 * this.camera.zoom;
        const scaleY = height * 0.8 * this.camera.zoom;
        const wx = 0.5 - this.camera.x + (sx - width / 2) / scaleX;
        const wy = 0.5 - this.camera.y - (sy - height / 2) / scaleY;
        return { x: wx, y: wy };
    }

    drawGrid() {
        const { width, height } = this.canvas;
        this.ctx.strokeStyle = '#1e293b';
        this.ctx.setLineDash([2, 4]);
        this.ctx.lineWidth = 1;
        
        // Vertical grid lines & X-axis ticks
        this.ctx.fillStyle = '#94a3b8';
        this.ctx.font = '10px Inter, sans-serif';
        this.ctx.textAlign = 'center';
        
        for (let x = 0; x <= 1.0; x += 0.2) {
            const p1 = this.worldToScreen(x, 0.5);
            this.ctx.beginPath(); 
            this.ctx.moveTo(p1.x, 0); 
            this.ctx.lineTo(p1.x, height); 
            this.ctx.stroke();
            
            // X-axis label (parametric xi)
            this.ctx.fillText(x.toFixed(1), p1.x, height - 10);
        }
        this.ctx.setLineDash([]);

        // Horizontal middle axis (Y=0.5 in world, which is 0 displacement)
        const midline = this.worldToScreen(0, 0.5);
        this.ctx.strokeStyle = '#334155';
        this.ctx.lineWidth = 2;
        this.ctx.beginPath(); 
        this.ctx.moveTo(0, midline.y); 
        this.ctx.lineTo(width, midline.y); 
        this.ctx.stroke();

        // Axis Titles
        this.ctx.fillStyle = '#cbd5e1';
        this.ctx.font = 'bold 12px Inter, sans-serif';
        this.ctx.fillText('Parametric Coordinate (ξ)', width / 2, height - 25);
        
        this.ctx.save();
        this.ctx.translate(20, height / 2);
        this.ctx.rotate(-Math.PI / 2);
        this.ctx.fillText('Displacement / Basis Value', 0, 0);
        this.ctx.restore();
    }

    drawBasis(nurbs) {
        const colors = ['#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6'];
        const numBasis = nurbs.controlPoints.length;
        this.ctx.lineWidth = 2;
        for (let i = 0; i < numBasis; i++) {
            this.ctx.strokeStyle = colors[i % colors.length];
            this.ctx.beginPath();
            for (let step = 0; step <= 100; step++) {
                const xi = step / 100;
                const val = nurbs.evaluateAllBasis(xi)[i];
                const pt = this.worldToScreen(xi, 0.5 + val * 0.4); // Offset so basis is visible
                step === 0 ? this.ctx.moveTo(pt.x, pt.y) : this.ctx.lineTo(pt.x, pt.y);
            }
            this.ctx.stroke();
        }
    }

    drawCurve(nurbs, color = '#f8fafc') {
        this.ctx.strokeStyle = color;
        this.ctx.lineWidth = 5;
        this.ctx.beginPath();
        const steps = 200;
        for (let step = 0; step <= steps; step++) {
            const xi = step / steps;
            const ptWorld = nurbs.evaluate(xi);
            const ptScreen = this.worldToScreen(ptWorld.x, ptWorld.y);
            step === 0 ? this.ctx.moveTo(ptScreen.x, ptScreen.y) : this.ctx.lineTo(ptScreen.x, ptScreen.y);
        }
        this.ctx.stroke();
    }

    drawStressGradient(nurbs, physicsState) {
        const { stresses, maxStress } = physicsState;
        const steps = stresses.length - 1;
        this.ctx.lineWidth = 6;
        
        for (let i = 0; i < steps; i++) {
            const xi1 = i / steps;
            const xi2 = (i + 1) / steps;
            const s1 = stresses[i];
            const s2 = stresses[i+1];
            const avgStress = (s1 + s2) / 2;
            
            // Map stress to color (Blue-White-Red)
            const normalized = maxStress > 0 ? avgStress / maxStress : 0; // -1 to 1
            this.ctx.strokeStyle = this.getStressColor(normalized);
            
            const pt1World = nurbs.evaluate(xi1);
            const pt2World = nurbs.evaluate(xi2);
            const pt1 = this.worldToScreen(pt1World.x, pt1World.y);
            const pt2 = this.worldToScreen(pt2World.x, pt2World.y);
            
            this.ctx.beginPath();
            this.ctx.moveTo(pt1.x, pt1.y);
            this.ctx.lineTo(pt2.x, pt2.y);
            this.ctx.stroke();
        }
    }

    getStressColor(val) {
        // val is -1 to 1. -1: Blue, 0: White/Gray, 1: Red
        const r = val > 0 ? 244 : Math.floor(255 + val * 200);
        const g = val > 0 ? Math.floor(255 - val * 200) : Math.floor(255 + val * 200);
        const b = val < 0 ? 255 : Math.floor(255 - val * 200);
        return `rgb(${r}, ${g}, ${b})`;
    }

    drawReferenceFEM(femData, visualScale, physicsMode = 'bending') {
        this.ctx.strokeStyle = '#10b981';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([5, 5]);
        this.ctx.beginPath();
        femData.forEach((pt, i) => {
            let wx = pt.x, wy = 0.5 - pt.y * visualScale;
            if (physicsMode === 'axial') {
                wx = pt.x + pt.y * visualScale;
                wy = 0.5;
            }
            const pScreen = this.worldToScreen(wx, wy);
            i === 0 ? this.ctx.moveTo(pScreen.x, pScreen.y) : this.ctx.lineTo(pScreen.x, pScreen.y);
        });
        this.ctx.stroke();
        this.ctx.setLineDash([]);
        
        // Nodes
        this.ctx.fillStyle = 'rgba(16, 185, 129, 0.5)';
        femData.forEach((pt, i) => {
            if (i % 10 === 0) {
                let wx = pt.x, wy = 0.5 - pt.y * visualScale;
                if (physicsMode === 'axial') {
                    wx = pt.x + pt.y * visualScale;
                    wy = 0.5;
                }
                const pScreen = this.worldToScreen(wx, wy);
                this.ctx.beginPath(); this.ctx.arc(pScreen.x, pScreen.y, 3, 0, Math.PI*2); this.ctx.fill();
            }
        });
    }

    drawAnalyticalCurve(params) {
        const { loadPos, effectiveLoad, engine, visualScale } = params;
        this.ctx.strokeStyle = '#facc15';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([2, 4]);
        this.ctx.beginPath();
        for (let step = 0; step <= 100; step++) {
            const xi = step / 100;
            const defl = engine.getAnalyticalBeamDeflection(xi, loadPos, effectiveLoad);
            const wy = 0.5 - defl * visualScale;
            const pt = this.worldToScreen(xi, wy);
            step === 0 ? this.ctx.moveTo(pt.x, pt.y) : this.ctx.lineTo(pt.x, pt.y);
        }
        this.ctx.stroke();
        this.ctx.setLineDash([]);
    }

    drawCrosshair() {
        const { x, y } = this.mousePos;
        if (x < 0 || y < 0 || x > this.canvas.width || y > this.canvas.height) return;
        this.ctx.strokeStyle = 'rgba(148, 163, 184, 0.2)';
        this.ctx.beginPath();
        this.ctx.moveTo(x, 0); this.ctx.lineTo(x, this.canvas.height);
        this.ctx.moveTo(0, y); this.ctx.lineTo(this.canvas.width, y);
        this.ctx.stroke();
    }
}
