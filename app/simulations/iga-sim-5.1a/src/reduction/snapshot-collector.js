/**
 * Phase 5.1a - Snapshot Collector
 * Manages the storage and retrieval of displacement fields during the sweep.
 */

class SnapshotCollector {
    constructor() {
        this.snapshots = [];
        this.parameterLog = [];
    }

    add(params, displacement) {
        // Store displacement field (Float64Array)
        this.snapshots.push(new Float64Array(displacement));
        
        // Log parameters
        this.parameterLog.push({ ...params });
        
        console.log(`[Collector] Added snapshot for R=${params.R}, L1=${params.L1}, L2=${params.L2}`);
    }

    clear() {
        this.snapshots = [];
        this.parameterLog = [];
    }

    get count() {
        return this.snapshots.length;
    }
}

window.SnapshotCollector = SnapshotCollector;
