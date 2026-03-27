import { BasisFunctions, NURBS, Curve, createCircle, hRefine, pRefine, kRefine } from './engine/iga.js';

class CircleApp {
    constructor() {
        this.canvas = document.getElementById('sim-canvas');
        this.ctx = this.canvas.getContext('2d');

        // Initialize with standard circle
        this.original = createCircle(0, 0, 1);
        this.current = { ...this.original, points: [...this.original.points], weights: [...this.original.weights], knots: [...this.original.knots] };

        this.refinementHistory = [];

        this.showControlPolygon = true;
        this.showBasis = false;
        this.showKnots = true;

        this.setupEventListeners();
        this.resize();
        window.addEventListener('resize', () => this.resize());

        this.updateStats();
        this.render();
    }

    resize() {
        this.canvas.width = this.canvas.parentElement.clientWidth;
        this.canvas.height = this.canvas.parentElement.clientHeight;
        this.render();
    }

    worldToScreen(x, y) {
        const cx = this.canvas.width / 2;
        const cy = this.canvas.height / 2;
        const scale = Math.min(cx, cy) * 0.35;
        return { x: cx + x * scale, y: cy - y * scale };
    }

    screenToWorld(sx, sy) {
        const cx = this.canvas.width / 2;
        const cy = this.canvas.height / 2;
        const scale = Math.min(cx, cy) * 0.35;
        return { x: (sx - cx) / scale, y: -(sy - cy) / scale };
    }

    applyRefinement(type) {
        let { degree, knots, points, weights } = this.current;
        let result;

        switch (type) {
            case 'p': {
                const elevateBy = parseInt(document.getElementById('slider-p-elevate').value);
                result = { degree, knots, points, weights };
                for (let i = 0; i < elevateBy; i++) {
                    result = pRefine(result.degree, result.knots, result.points, result.weights);
                }
                this.refinementHistory.push(`P+${elevateBy}`);
                break;
            }
            case 'h': {
                const insertions = parseInt(document.getElementById('slider-h-insert').value);
                result = { degree, knots, points, weights };
                for (let i = 0; i < insertions; i++) {
                    result = hRefine(result.degree, result.knots, result.points, result.weights);
                }
                this.refinementHistory.push(`H×${insertions}`);
                break;
            }
            case 'k': {
                const param = parseInt(document.getElementById('slider-k-param').value);
                // k-refinement: elevate by param, then refine once
                result = { degree, knots, points, weights };
                for (let i = 0; i < param; i++) {
                    result = pRefine(result.degree, result.knots, result.points, result.weights);
                }
                result = hRefine(result.degree, result.knots, result.points, result.weights);
                this.refinementHistory.push(`K(+${param},h)`);
                break;
            }
        }

        if (result) {
            this.current = result;
            this.updateStats();
            this.updateHistoryUI();
            this.render();
        }
    }

    reset() {
        this.current = {
            ...this.original,
            points: [...this.original.points],
            weights: [...this.original.weights],
            knots: [...this.original.knots]
        };
        this.refinementHistory = [];
        this.updateStats();
        this.updateHistoryUI();
        this.render();
    }

    updateStats() {
        const { degree, points, knots } = this.current;
        const uniqueKnots = [...new Set(knots)];
        const numElements = uniqueKnots.length - 1;

        // Determine continuity: max interior multiplicity
        let maxMult = 0;
        for (let i = 1; i < uniqueKnots.length - 1; i++) {
            const mult = knots.filter(k => Math.abs(k - uniqueKnots[i]) < 1e-10).length;
            maxMult = Math.max(maxMult, mult);
        }
        const continuity = maxMult > 0 ? degree - maxMult : degree;

        document.getElementById('stat-degree').textContent = degree;
        document.getElementById('stat-cp').textContent = points.length;
        document.getElementById('stat-elements').textContent = numElements;
        document.getElementById('stat-continuity').textContent = `C${continuity}`;
        document.getElementById('stat-knots').textContent = knots.length;
    }

    updateHistoryUI() {
        const container = document.getElementById('history-list');
        if (this.refinementHistory.length === 0) {
            container.innerHTML = '<span style="color: var(--text-secondary); font-size: 0.8rem;">No refinements applied</span>';
            return;
        }
        container.innerHTML = this.refinementHistory.map((type, i) => {
            const colors = { P: '#a78bfa', H: '#34d399', K: '#f59e0b' };
            return `<span class="history-badge" style="background: ${colors[type]}20; color: ${colors[type]}; border: 1px solid ${colors[type]}40;">${i + 1}. ${type}-refine</span>`;
        }).join('');
    }

    drawGrid() {
        const ctx = this.ctx;
        ctx.strokeStyle = '#1e293b';
        ctx.lineWidth = 1;

        // Draw grid lines
        for (let i = -3; i <= 3; i++) {
            const from = this.worldToScreen(i, -3);
            const to = this.worldToScreen(i, 3);
            ctx.beginPath();
            ctx.moveTo(from.x, from.y);
            ctx.lineTo(to.x, to.y);
            ctx.stroke();

            const from2 = this.worldToScreen(-3, i);
            const to2 = this.worldToScreen(3, i);
            ctx.beginPath();
            ctx.moveTo(from2.x, from2.y);
            ctx.lineTo(to2.x, to2.y);
            ctx.stroke();
        }

        // Draw axes
        ctx.strokeStyle = '#334155';
        ctx.lineWidth = 2;
        const ox = this.worldToScreen(-3.2, 0), oxe = this.worldToScreen(3.2, 0);
        ctx.beginPath(); ctx.moveTo(ox.x, ox.y); ctx.lineTo(oxe.x, oxe.y); ctx.stroke();
        const oy = this.worldToScreen(0, -3.2), oye = this.worldToScreen(0, 3.2);
        ctx.beginPath(); ctx.moveTo(oy.x, oy.y); ctx.lineTo(oye.x, oye.y); ctx.stroke();

        // Axis Labels
        ctx.font = 'bold 14px "Inter"';
        ctx.fillStyle = '#94a3b8';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        
        const xPos = this.worldToScreen(3.4, 0);
        ctx.fillText('X', xPos.x, xPos.y);
        
        const yPos = this.worldToScreen(0, 3.4);
        ctx.fillText('Y', yPos.x, yPos.y);

        // Numeric Ticks
        ctx.font = '10px "JetBrains Mono"';
        ctx.fillStyle = '#64748b';
        for (let i = -3; i <= 3; i++) {
            if (i === 0) continue;
            // X-Ticks
            const tx = this.worldToScreen(i, 0);
            ctx.beginPath(); ctx.moveTo(tx.x, tx.y - 4); ctx.lineTo(tx.x, tx.y + 4); ctx.stroke();
            ctx.fillText(i.toString(), tx.x, tx.y + 15);

            // Y-Ticks
            const ty = this.worldToScreen(0, i);
            ctx.beginPath(); ctx.moveTo(ty.x - 4, ty.y); ctx.lineTo(ty.x + 4, ty.y); ctx.stroke();
            ctx.textAlign = 'right';
            ctx.fillText(i.toString(), ty.x - 10, ty.y);
            ctx.textAlign = 'center';
        }
    }

    drawControlPolygon() {
        if (!this.showControlPolygon) return;
        const ctx = this.ctx;
        const pts = this.current.points;

        // Lines
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.3)';
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        for (let i = 0; i < pts.length; i++) {
            const s = this.worldToScreen(pts[i].x, pts[i].y);
            if (i === 0) ctx.moveTo(s.x, s.y);
            else ctx.lineTo(s.x, s.y);
        }
        ctx.stroke();
        ctx.setLineDash([]);

        // Points
        for (let i = 0; i < pts.length; i++) {
            const s = this.worldToScreen(pts[i].x, pts[i].y);
            const w = this.current.weights[i];
            const isWeighted = Math.abs(w - 1.0) > 0.01;

            ctx.beginPath();
            ctx.arc(s.x, s.y, isWeighted ? 6 : 5, 0, Math.PI * 2);
            ctx.fillStyle = isWeighted ? '#f59e0b' : '#10b981';
            ctx.fill();
            ctx.strokeStyle = 'rgba(255,255,255,0.4)';
            ctx.lineWidth = 1.5;
            ctx.stroke();

            // Weight label for weighted CPs
            if (isWeighted) {
                ctx.font = '10px "JetBrains Mono", monospace';
                ctx.fillStyle = '#f59e0b';
                ctx.textAlign = 'center';
                ctx.fillText(`w=${w.toFixed(3)}`, s.x, s.y - 12);
            }
        }
    }

    drawCurve() {
        const ctx = this.ctx;
        const { degree, knots, points, weights } = this.current;
        const steps = 300;

        ctx.beginPath();
        for (let i = 0; i <= steps; i++) {
            const xi = i / steps;
            const pt = Curve.evaluate(xi, degree, knots, points, weights);
            const s = this.worldToScreen(pt.x, pt.y);
            if (i === 0) ctx.moveTo(s.x, s.y);
            else ctx.lineTo(s.x, s.y);
        }
        ctx.strokeStyle = '#3b82f6';
        ctx.lineWidth = 2.5;
        ctx.shadowColor = 'rgba(59, 130, 246, 0.5)';
        ctx.shadowBlur = 10;
        ctx.stroke();
        ctx.shadowBlur = 0;

        // Draw knot positions on curve
        if (this.showKnots) {
            const uniqueKnots = [...new Set(knots)];
            for (const u of uniqueKnots) {
                const pt = Curve.evaluate(Math.min(u, 0.9999), degree, knots, points, weights);
                const s = this.worldToScreen(pt.x, pt.y);
                const mult = knots.filter(k => Math.abs(k - u) < 1e-10).length;

                ctx.beginPath();
                ctx.arc(s.x, s.y, 4, 0, Math.PI * 2);
                ctx.fillStyle = '#ef4444';
                ctx.fill();

                if (mult > 1) {
                    ctx.font = '9px "JetBrains Mono", monospace';
                    ctx.fillStyle = '#ef4444';
                    ctx.textAlign = 'center';
                    ctx.fillText(`m=${mult}`, s.x, s.y + 16);
                }
            }
        }
    }

    drawBasisFunctions() {
        if (!this.showBasis) return;

        const ctx = this.ctx;
        const { degree, knots, points, weights } = this.current;
        const n = points.length;
        const steps = 200;

        const colors = [
            '#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6',
            '#ec4899', '#06b6d4', '#84cc16', '#f97316', '#6366f1',
            '#14b8a6', '#e11d48', '#a3e635', '#0ea5e9', '#d946ef'
        ];

        // Draw in a strip at the bottom of the canvas
        const basisH = this.canvas.height * 0.25;
        const basisY = this.canvas.height - basisH - 20;
        const basisW = this.canvas.width - 80;
        const basisX = 40;

        // Background
        ctx.fillStyle = 'rgba(15, 23, 42, 0.85)';
        ctx.fillRect(basisX - 10, basisY - 10, basisW + 20, basisH + 30);
        ctx.strokeStyle = '#334155';
        ctx.lineWidth = 1;
        ctx.strokeRect(basisX - 10, basisY - 10, basisW + 20, basisH + 30);

        // Label
        ctx.font = '11px "Inter", sans-serif';
        ctx.fillStyle = '#94a3b8';
        ctx.textAlign = 'left';
        ctx.fillText('Basis Functions', basisX, basisY - 15);

        for (let i = 0; i < n; i++) {
            ctx.beginPath();
            ctx.strokeStyle = colors[i % colors.length];
            ctx.lineWidth = 1.5;

            for (let s = 0; s <= steps; s++) {
                const xi = s / steps;
                const N = BasisFunctions.evaluate(i, degree, xi, knots);
                const px = basisX + xi * basisW;
                const py = basisY + basisH - N * basisH;

                if (s === 0) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            }
            ctx.stroke();
        }
    }

    render() {
        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        this.drawGrid();
        this.drawCurve();
        this.drawControlPolygon();
        this.drawBasisFunctions();
    }

    setupEventListeners() {
        // Refinement Apply Buttons
        document.getElementById('btn-p-apply').addEventListener('click', () => this.applyRefinement('p'));
        document.getElementById('btn-h-apply').addEventListener('click', () => this.applyRefinement('h'));
        document.getElementById('btn-k-apply').addEventListener('click', () => this.applyRefinement('k'));
        document.getElementById('btn-reset').addEventListener('click', () => this.reset());

        // Sliders Real-time feedback
        document.getElementById('slider-p-elevate').addEventListener('input', (e) => {
            document.getElementById('val-p-elevate').textContent = `+${e.target.value}`;
        });
        document.getElementById('slider-h-insert').addEventListener('input', (e) => {
            document.getElementById('val-h-insert').textContent = e.target.value;
        });
        document.getElementById('slider-k-param').addEventListener('input', (e) => {
            document.getElementById('val-k-param').textContent = `+${e.target.value}, ${e.target.value}`;
        });

        // Toggles
        document.getElementById('toggle-cp').addEventListener('change', (e) => {
            this.showControlPolygon = e.target.checked;
            this.render();
        });

        document.getElementById('toggle-basis').addEventListener('change', (e) => {
            this.showBasis = e.target.checked;
            this.render();
        });

        document.getElementById('toggle-knots').addEventListener('change', (e) => {
            this.showKnots = e.target.checked;
            this.render();
        });
    }
}

// Start
window.addEventListener('DOMContentLoaded', () => new CircleApp());
