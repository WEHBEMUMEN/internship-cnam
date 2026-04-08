/**
 * Performance Monitoring Utility
 * Tracks FPS and CPU Processing Time
 */
class PerformanceMonitor {
    constructor() {
        this.fps = 0;
        this.cpuTime = 0;
        this.lastTime = performance.now();
        this.frameCount = 0;
        this.startTime = 0;
        
        this.createOverlay();
        this.tick();
    }

    createOverlay() {
        const overlay = document.createElement('div');
        overlay.id = 'performance-overlay';
        overlay.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            background: rgba(15, 23, 42, 0.8);
            backdrop-filter: blur(8px);
            border: 1px solid rgba(59, 130, 246, 0.3);
            border-radius: 12px;
            padding: 12px 16px;
            color: #f8fafc;
            font-family: 'JetBrains Mono', 'Monaco', monospace;
            font-size: 11px;
            z-index: 9999;
            pointer-events: none;
            display: flex;
            flex-direction: column;
            gap: 4px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.3);
        `;
        
        overlay.innerHTML = `
            <div style="display: flex; justify-content: space-between; gap: 20px;">
                <span style="color: #94a3b8;">FPS</span>
                <span id="perf-fps" style="color: #10b981; font-weight: 700;">0</span>
            </div>
            <div style="display: flex; justify-content: space-between; gap: 20px;">
                <span style="color: #94a3b8;">CPU</span>
                <span id="perf-cpu" style="color: #3b82f6; font-weight: 700;">0.00 ms</span>
            </div>
        `;
        
        document.body.appendChild(overlay);
    }

    startMeasure() {
        this.startTime = performance.now();
    }

    endMeasure() {
        const duration = performance.now() - this.startTime;
        this.cpuTime = duration;
        document.getElementById('perf-cpu').textContent = `${duration.toFixed(2)} ms`;
    }

    tick() {
        const now = performance.now();
        this.frameCount++;
        
        if (now - this.lastTime >= 1000) {
            this.fps = Math.round((this.frameCount * 1000) / (now - this.lastTime));
            document.getElementById('perf-fps').textContent = this.fps;
            this.frameCount = 0;
            this.lastTime = now;
        }
        
        requestAnimationFrame(() => this.tick());
    }
}

// Initialize on load
window.addEventListener('DOMContentLoaded', () => {
    window.perfMonitor = new PerformanceMonitor();
});
