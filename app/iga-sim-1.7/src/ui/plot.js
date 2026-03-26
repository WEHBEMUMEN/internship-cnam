/**
 * IGA Plot Renderer
 * High-performance Canvas rendering for basis functions
 */

export class BasisPlot {
    constructor(canvasOrId) {
        this.canvas = typeof canvasOrId === 'string' ? document.getElementById(canvasOrId) : canvasOrId;
        this.ctx = this.canvas.getContext('2d');
        
        // Camera state
        this.camera = { x: 0, y: 0, zoom: 1.0 };
        this.minZoom = 0.5;
        this.maxZoom = 20.0;
        
        // Interaction state
        this.isDragging = false;
        this.lastMouse = { x: 0, y: 0 };
        this.mousePos = { x: -1000, y: -1000, worldX: 0, worldY: 0 }; // far away initially
        this.hoverValue = null;
        
        this.resize();
        
        // Use ResizeObserver handles manual `resize: both` handles AND window resizing automatically
        const wrapper = this.canvas.parentElement;
        this.resizeObserver = new ResizeObserver(() => this.resize());
        this.resizeObserver.observe(wrapper);
        
        this.setupEvents();
    }

    setupEvents() {
        this.canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            this.handleZoom(e.deltaY < 0 ? 1.1 : 0.9, e.offsetX, e.offsetY);
        });

        this.canvas.addEventListener('mousedown', (e) => {
            // Pan with right click or middle click
            if (e.button === 2 || e.button === 1 || e.button === 0) {
                this.isDragging = true;
                this.lastMouse = { x: e.offsetX, y: e.offsetY };
                this.canvas.style.cursor = 'grabbing';
            }
        });

        window.addEventListener('mouseup', () => {
            this.isDragging = false;
            this.canvas.style.cursor = 'crosshair';
        });

        this.canvas.addEventListener('mousemove', (e) => {
            const rect = this.canvas.getBoundingClientRect();
            this.mousePos.x = e.clientX - rect.left;
            this.mousePos.y = e.clientY - rect.top;
            
            const wPos = this.screenToWorld(this.mousePos.x, this.mousePos.y);
            this.mousePos.worldX = wPos.x;
            this.mousePos.worldY = wPos.y;

            if (this.isDragging) {
                const dx = e.offsetX - this.lastMouse.x;
                const dy = e.offsetY - this.lastMouse.y;
                
                const { width, height } = this.canvas;
                const cw = width * 0.8 * this.camera.zoom;
                const ch = height * 0.8 * this.camera.zoom;
                
                this.camera.x -= dx / cw;
                this.camera.y += dy / ch; // + because screen Y is inverted
                
                this.lastMouse = { x: e.offsetX, y: e.offsetY };
            }
            
            // Update UI readout if present
            const readout = document.getElementById('coord-readout');
            if (readout) {
                readout.textContent = `ξ: ${this.mousePos.worldX.toFixed(3)} | y: ${this.mousePos.worldY.toFixed(3)}`;
            }
        });
        
        this.canvas.addEventListener('contextmenu', e => e.preventDefault());
    }

    handleZoom(factor, mx = this.canvas.width/2, my = this.canvas.height/2) {
        const oldZoom = this.camera.zoom;
        let newZoom = oldZoom * factor;
        newZoom = Math.max(this.minZoom, Math.min(newZoom, this.maxZoom));
        
        // Center zoom on mouse cursor
        const { width, height } = this.canvas;
        const cw = width * 0.8;
        const ch = height * 0.8;
        
        const dx = mx - width / 2;
        const dy = height / 2 - my;
        
        // The mouse position in world should remain constant before and after zoom
        const wx = this.camera.x + dx / (cw * oldZoom);
        const wy = this.camera.y + dy / (ch * oldZoom);
        
        this.camera.zoom = newZoom;
        
        this.camera.x = wx - dx / (cw * newZoom);
        this.camera.y = wy - dy / (ch * newZoom);
    }

    resetView() {
        this.camera = { x: 0, y: 0, zoom: 1.0 };
    }

    worldToScreen(wx, wy) {
        const { width, height } = this.canvas;
        const cw = width * 0.8;
        const ch = height * 0.8;
        const cx = this.camera.x + 0.5;
        const cy = this.camera.y + 0.5;

        const sx = width / 2 + (wx - cx) * cw * this.camera.zoom;
        const sy = height / 2 - (wy - cy) * ch * this.camera.zoom;

        return { x: sx, y: sy };
    }

    screenToWorld(sx, sy) {
        const { width, height } = this.canvas;
        const cw = width * 0.8;
        const ch = height * 0.8;
        const cx = this.camera.x + 0.5;
        const cy = this.camera.y + 0.5;

        const wx = cx + (sx - width / 2) / (cw * this.camera.zoom);
        const wy = cy - (sy - height / 2) / (ch * this.camera.zoom);

        return { x: wx, y: wy };
    }

    resize() {
        const rect = this.canvas.parentElement.getBoundingClientRect();
        this.canvas.width = rect.width;
        this.canvas.height = rect.height;
    }

    clear() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    }

    drawGrid() {
        const { width, height } = this.canvas;
        this.ctx.lineWidth = 1;
        this.ctx.font = '12px "JetBrains Mono", monospace';
        this.ctx.textAlign = 'center';
        
        // Adaptive grid step
        let step = 0.25;
        if (this.camera.zoom > 3) step = 0.1;
        if (this.camera.zoom > 8) step = 0.05;
        
        // Find screen extents in world
        const tr = this.screenToWorld(width, 0);
        const bl = this.screenToWorld(0, height);
        const wMinX = Math.floor(bl.x / step) * step;
        const wMaxX = Math.ceil(tr.x / step) * step;
        const wMinY = Math.floor(bl.y / step) * step;
        const wMaxY = Math.ceil(tr.y / step) * step;

        this.ctx.strokeStyle = '#334155';
        this.ctx.fillStyle = '#64748b';

        // X Grid
        for (let x = wMinX; x <= wMaxX; x += step) {
            const p1 = this.worldToScreen(x, bl.y);
            const p2 = this.worldToScreen(x, tr.y);
            this.ctx.beginPath(); this.ctx.moveTo(p1.x, 0); this.ctx.lineTo(p2.x, height); this.ctx.stroke();
            
            // X Axis ticks
            const pX = this.worldToScreen(x, 0);
            if (pX.x > 0 && pX.x < width) {
                this.ctx.fillText(x.toFixed(step >= 0.1 ? 1 : 2), pX.x, pX.y + 15);
            }
        }

        // Y Grid
        for (let y = wMinY; y <= wMaxY; y += step) {
            const p1 = this.worldToScreen(bl.x, y);
            const p2 = this.worldToScreen(tr.x, y);
            this.ctx.beginPath(); this.ctx.moveTo(0, p1.y); this.ctx.lineTo(width, p2.y); this.ctx.stroke();
            
            // Y Axis ticks
            const pY = this.worldToScreen(0, y);
            if (pY.y > 0 && pY.y < height && Math.abs(y) > 0.01) {
                this.ctx.textAlign = 'right';
                this.ctx.fillText(y.toFixed(step >= 0.1 ? 1 : 2), pY.x - 8, pY.y + 4);
                this.ctx.textAlign = 'center';
            }
        }

        // Highlight Axes
        this.ctx.strokeStyle = '#475569';
        this.ctx.lineWidth = 2;
        const o1 = this.worldToScreen(0, bl.y); const o2 = this.worldToScreen(0, tr.y);
        this.ctx.beginPath(); this.ctx.moveTo(o1.x, 0); this.ctx.lineTo(o2.x, height); this.ctx.stroke(); // Y Axis
        
        const x1 = this.worldToScreen(bl.x, 0); const x2 = this.worldToScreen(tr.x, 0);
        this.ctx.beginPath(); this.ctx.moveTo(0, x1.y); this.ctx.lineTo(width, x2.y); this.ctx.stroke(); // X Axis
    }

    drawCrosshair() {
        const { width, height } = this.canvas;
        const { x, y } = this.mousePos;
        
        if (x < 0 || y < 0 || x > width || y > height) return;

        this.ctx.strokeStyle = '#94a3b8';
        this.ctx.lineWidth = 1;
        this.ctx.setLineDash([4, 4]);

        this.ctx.beginPath();
        this.ctx.moveTo(x, 0);
        this.ctx.lineTo(x, height);
        this.ctx.moveTo(0, y);
        this.ctx.lineTo(width, y);
        this.ctx.stroke();
        this.ctx.setLineDash([]);
        
        // Tooltip drawing
        if (this.hoverValue !== null) {
            this.ctx.fillStyle = '#0f172a';
            this.ctx.strokeStyle = '#3b82f6';
            this.ctx.lineWidth = 1;
            
            const txt = `def: ${this.hoverValue.toFixed(4)}`;
            const tw = this.ctx.measureText(txt).width;
            
            this.ctx.beginPath();
            this.ctx.roundRect(x + 10, y - 25, tw + 16, 24, 4);
            this.ctx.fill();
            this.ctx.stroke();
            
            this.ctx.fillStyle = '#f8fafc';
            this.ctx.textAlign = 'left';
            this.ctx.fillText(txt, x + 18, y - 9);
            this.ctx.textAlign = 'center';
        }
    }

    drawBasis(nurbs) {
        const colors = ['#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#06b6d4', '#f97316', '#22c55e', '#a855f7'];
        const numBasis = nurbs.controlPoints.length;
        
        for (let i = 0; i < numBasis; i++) {
            const color = colors[i % colors.length];
            this.ctx.strokeStyle = color;
            this.ctx.lineWidth = 3;
            this.ctx.beginPath();

            let first = true;
            const segments = 100;
            for (let step = 0; step <= segments; step++) {
                const xi = step / segments;
                const val = nurbs.evaluateAllBasis(xi)[i];
                const pt = this.worldToScreen(xi, val);
                if (first) { this.ctx.moveTo(pt.x, pt.y); first = false; }
                else { this.ctx.lineTo(pt.x, pt.y); }
            }
            this.ctx.stroke();
            
            // Fill area under curve
            this.ctx.fillStyle = color + '22';
            const base1 = this.worldToScreen(1, 0);
            const base0 = this.worldToScreen(0, 0);
            this.ctx.lineTo(base1.x, base1.y);
            this.ctx.lineTo(base0.x, base0.y);
            this.ctx.fill();
        }
    }

    /**
     * Draws the physical NURBS curve C(xi) - The ROM Solution
     */
    drawCurve(nurbs, color = '#3b82f6', visualScale = 0.05) {
        this.ctx.strokeStyle = color;
        this.ctx.lineWidth = 4;
        this.ctx.lineJoin = 'round';
        this.ctx.setLineDash([]);

        this.ctx.beginPath();
        let first = true;

        const segments = 200;
        for (let step = 0; step <= segments; step++) {
            const xi = step / segments;
            const ptWorld = nurbs.evaluate(xi);
            
            // Map physical deflection to visualization space
            const wy = 0.5 - ptWorld.y * visualScale;
            const ptScreen = this.worldToScreen(xi, wy);

            if (first) { this.ctx.moveTo(ptScreen.x, ptScreen.y); first = false; }
            else { this.ctx.lineTo(ptScreen.x, ptScreen.y); }
        }
        this.ctx.stroke();
    }

    drawControlPoints(nurbs) {
        nurbs.controlPoints.forEach((pt) => {
            const pScreen = this.worldToScreen(pt.x, pt.y);
            this.ctx.fillStyle = '#94a3b8';
            this.ctx.beginPath();
            this.ctx.arc(pScreen.x, pScreen.y, 5, 0, Math.PI * 2);
            this.ctx.fill();
        });
    }

    /**
     * Draws the High-Resolution Reference FEM Solution (Ground Truth)
     */
    drawReferenceFEM(femData, visualScale = 0.05) {
        this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([5, 5]);

        this.ctx.beginPath();
        let first = true;

        femData.forEach(pt => {
            const wy = 0.5 - pt.y * visualScale;
            const ptScreen = this.worldToScreen(pt.x, wy);
            if (first) {
                this.ctx.moveTo(ptScreen.x, ptScreen.y);
                first = false;
            } else {
                this.ctx.lineTo(ptScreen.x, ptScreen.y);
            }
        });
        this.ctx.stroke();
        this.ctx.setLineDash([]);

        // Highlight every 10th FEM node as a ground truth marker
        this.ctx.fillStyle = 'rgba(255, 255, 255, 0.3)';
        femData.forEach((pt, i) => {
            if (i % 10 === 0) {
                const wy = 0.5 - pt.y * visualScale;
                const ptScreen = this.worldToScreen(pt.x, wy);
                this.ctx.beginPath();
                this.ctx.arc(ptScreen.x, ptScreen.y, 3, 0, Math.PI * 2);
                this.ctx.fill();
            }
        });
    }

    drawAnalyticalCurve(physicsParams) {
        const { loadPos, effectiveLoad, engine, visualScale = 0.05 } = physicsParams;
        this.ctx.strokeStyle = '#facc15'; // yellow color
        this.ctx.lineWidth = 3;
        this.ctx.setLineDash([2, 6]);

        this.ctx.beginPath();
        let first = true;
        const segments = 200;

        for (let step = 0; step <= segments; step++) {
            const xi = step / segments;
            const defl = engine.getAnalyticalBeamDeflection(xi, loadPos, effectiveLoad);
            
            // Map true beam deflection to our normalized Y axis 0.5 center.
            const wy = 0.5 - defl * visualScale; 
            const ptScreen = this.worldToScreen(xi, wy);

            if (first) { this.ctx.moveTo(ptScreen.x, ptScreen.y); first = false; }
            else { this.ctx.lineTo(ptScreen.x, ptScreen.y); }
        }
        
        this.ctx.stroke();
        this.ctx.setLineDash([]);
    }
}
