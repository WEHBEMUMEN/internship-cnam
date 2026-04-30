/**
 * ECSW Audit & Diagnostic Tool
 * Generates technical reports to verify hyper-reduction accuracy.
 */
window.ECSWAudit = {
    /**
     * Generate a comprehensive console report for ECSW training.
     */
    report: function(trainedData, pods) {
        const { sampleElements, k, nDofs, weights } = trainedData;
        const totalElements = weights.length;
        const activeCount = sampleElements.length;
        const sparsity = ((1 - activeCount / totalElements) * 100).toFixed(1);

        console.log(`%c════════════════════════════════════════════════════════`, 'color: #8b5cf6; font-weight: bold');
        console.log(`%c  ECSW ENGINE AUDIT REPORT`, 'color: #8b5cf6; font-weight: bold; font-size: 1.1rem');
        console.log(`%c  Generated: ${new Date().toISOString()}`, 'color: #64748b; font-style: italic');
        console.log(`%c════════════════════════════════════════════════════════`, 'color: #8b5cf6; font-weight: bold');

        console.log(`%c0. CONFIGURATION`, 'color: #1e293b; font-weight: bold');
        console.log(`   POD Modes (k)      : ${k}`);
        console.log(`   Total Elements     : ${totalElements}`);
        console.log(`   System DOFs        : ${nDofs}`);
        console.log(``);

        console.log(`%c1. SPARSITY & PERFORMANCE`, 'color: #1e293b; font-weight: bold');
        console.log(`   Active Elements    : ${activeCount}`);
        console.log(`   Reduction Ratio    : ${sparsity}%`);
        console.log(`   Speedup Potential  : ~${(totalElements / activeCount).toFixed(1)}x (Assembly phase)`);
        console.log(``);

        console.log(`%c2. WEIGHT DISTRIBUTION`, 'color: #1e293b; font-weight: bold');
        // Simple histogram logic
        const activeWeights = sampleElements.map(el => el.weight);
        const maxW = Math.max(...activeWeights);
        const minW = Math.min(...activeWeights);
        console.log(`   Weight Range       : [${minW.toFixed(4)}, ${maxW.toFixed(4)}]`);
        
        // Show top 5 weights
        const sorted = [...sampleElements].sort((a,b) => b.weight - a.weight).slice(0, 5);
        console.log(`   Top 5 Contributions:`);
        sorted.forEach((el, i) => {
            console.log(`     ${i+1}. Element[${el.i},${el.j}] : weight = ${el.weight.toFixed(4)}`);
        });
        console.log(``);

        console.log(`%c3. BASIS ENERGY (POD)`, 'color: #1e293b; font-weight: bold');
        if (pods && pods.singularValues) {
            const sv = pods.singularValues;
            const totalEnergy = sv.reduce((a, b) => a + b, 0);
            const captured = sv.slice(0, k).reduce((a, b) => a + b, 0);
            console.log(`   Captured Energy    : ${(captured / totalEnergy * 100).toFixed(6)}%`);
            console.log(`   Last SV (k=${k})   : ${sv[k-1].toExponential(4)}`);
        } else {
            console.log(`   [Basis data not provided to audit]`);
        }

        console.log(`%c════════════════════════════════════════════════════════`, 'color: #8b5cf6; font-weight: bold');
    },

    /**
     * Log comparative internal force data during Newton iterations.
     */
    logIteration: function(iter, norm, relError) {
        const color = relError < 0.05 ? '#10b981' : (relError < 0.2 ? '#f59e0b' : '#ef4444');
        console.log(
            `%c[ECSW Iter ${iter}] %cResidual: ${norm.toExponential(4)} %c| Rel. Force Error: ${(relError * 100).toFixed(2)}%`,
            'color: #8b5cf6; font-weight: bold',
            'color: #1e293b',
            `color: ${color}; font-weight: bold`
        );
    }
};
