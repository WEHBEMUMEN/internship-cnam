/**
 * ECSW Setup & Integration
 * Manages the global engine state and provides the interface for the app.
 */
class ECSWManager {
    constructor() {
        this.trainer = new window.ECSWTrainer();
        this.online = new window.ECSWOnlineEngine();
        this.trainedData = null;
    }

    /**
     * High-level training interface.
     */
    async train(fomSolver, romEngine, patch, snapshots, tolerance) {
        this.trainedData = await this.trainer.train({
            fomSolver,
            romEngine,
            patch,
            snapshots,
            tolerance
        });

        // Initialize the online engine with the results
        this.online.init(this.trainedData, romEngine.Phi);

        // Generate Audit Report
        if (window.ECSWAudit) {
            window.ECSWAudit.report(this.trainedData, romEngine.podInfo);
        }

        return {
            elementCount: this.trainedData.sampleElements.length,
            totalElements: this.trainer.allElements.length
        };
    }

    /**
     * Solve using the online engine.
     */
    solveReduced(fomSolver, patch, loads, options) {
        if (!this.trainedData) {
            throw new Error("ECSW: Cannot solve before training.");
        }
        return this.online.solveReduced(fomSolver, patch, loads, options);
    }

    /**
     * Accessor for active elements (used by UI for visualization).
     */
    get sampleElements() {
        return this.trainedData ? this.trainedData.sampleElements : [];
    }
}

// Replace the monolithic ECSWEngine with our new modular Manager
window.ECSWEngine = ECSWManager;
