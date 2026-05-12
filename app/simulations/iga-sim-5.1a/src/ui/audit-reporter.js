/**
 * Phase 5.1a - Audit Reporter
 * Handles the display of training analytics and system logs.
 */

class AuditReporter {
    constructor(logId) {
        this.logElement = document.getElementById(logId);
    }

    log(type, message) {
        const entry = document.createElement('div');
        entry.className = `log-entry ${type}`;
        
        const timestamp = new Date().toLocaleTimeString();
        let icon = '<i class="fa-solid fa-info-circle"></i>';
        
        if (type === 'system') icon = '<i class="fa-solid fa-microchip"></i>';
        if (type === 'success') icon = '<i class="fa-solid fa-check-circle" style="color: #10b981;"></i>';
        if (type === 'warning') icon = '<i class="fa-solid fa-triangle-exclamation" style="color: #f59e0b;"></i>';
        
        entry.innerHTML = `
            <span class="log-time">[${timestamp}]</span>
            <span class="log-icon">${icon}</span>
            <span class="log-msg">${message}</span>
        `;
        
        this.logElement.prepend(entry);
    }

    clear() {
        this.logElement.innerHTML = '';
    }

    reportStats(stats) {
        if (document.getElementById('stat-dofs')) {
            document.getElementById('stat-dofs').textContent = stats.nDofs;
        }
        if (document.getElementById('stat-snaps')) {
            document.getElementById('stat-snaps').textContent = stats.nSnaps;
        }
        
        // Reduction Ratio
        if (document.getElementById('stat-ratio')) {
            const ratio = ((stats.nDofs - stats.k) / stats.nDofs) * 100;
            document.getElementById('stat-ratio').textContent = `${ratio.toFixed(1)}%`;
        }

        // Error Handling
        const errorContainer = document.getElementById('stat-rom-container');
        if (errorContainer) {
            if (stats.mode === 'rom') {
                errorContainer.style.display = 'block';
                const errVal = document.getElementById('stat-error');
                if (errVal) errVal.textContent = `${(stats.error * 100).toFixed(4)}%`;
            } else {
                errorContainer.style.display = 'none';
            }
        }
    }
}

window.AuditReporter = AuditReporter;
