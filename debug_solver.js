const NURBS2D = require('./app/phase-2-core/src/nurbs-2d.js');
const { IGA2DSolver, GaussQuadrature2D } = require('./app/phase-2-core/src/iga-solver-2d.js');

const NURBSPresets = require('./app/phase-2-core/src/nurbs-presets.js');

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);
const patch = NURBSPresets.generateSheet();

const nU = patch.controlPoints.length;
const nV = patch.controlPoints[0].length;
const bcs = [];
for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });

const loads = [{ i: nU - 1, j: Math.floor(nV / 2), fx: 0, fy: -5000 }];

console.log("Running Linear...");
const uLin = solver.solve(patch, bcs, loads);
let maxLin = 0;
for(let i=0; i<uLin.length; i+=2) maxLin = Math.max(maxLin, Math.sqrt(uLin[i]**2 + uLin[i+1]**2));
console.log("Linear Max Deflection:", maxLin);

console.log("Running Nonlinear...");
const result = solver.solveNonlinear(patch, bcs, loads, { steps: 2, iterations: 10 });
console.log("Residual History:", result.residualHistory);
let maxNL = 0;
for(let i=0; i<result.u.length; i+=2) maxNL = Math.max(maxNL, Math.sqrt(result.u[i]**2 + result.u[i+1]**2));
console.log("Nonlinear Max Deflection:", maxNL);
