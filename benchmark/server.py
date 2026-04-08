from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import nutils, numpy
from nutils import mesh, function, solver
from nutils.expression_v2 import Namespace
from typing import List, Optional

app = FastAPI()

# Enable CORS so the web browser (e.g. from Netlify or localhost) can talk to this PC server
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SolveRequest(BaseModel):
    radius: float
    length: float
    E: float
    nu: float
    traction: float
    nrefine: int

@app.get("/")
def health():
    return {"status": "online", "engine": "Nutils IGA Core", "version": "1.0.0"}

@app.post("/solve")
async def solve_iga(req: SolveRequest):
    try:
        # 1. Topology
        # We use the standard 2-element base mesh (Mesh 1 in Fig 4.2)
        topo, (u, v) = mesh.rectilinear([numpy.linspace(0, 1, 3), numpy.linspace(0, 1, 2)])
        for _ in range(req.nrefine):
            topo = topo.refine(1)
            u, v = topo.sample('gauss', 1).eval([u, v]).T
            
        basis = topo.basis('spline', degree=2)
        
        # 2. Geometry (Hughes/Cottrell Benchmarking Mapping)
        t = numpy.tan(numpy.pi/8)
        w_arc = numpy.cos(numpy.pi/8)
        R = req.radius
        L = req.length
        
        def mid(p1, p2): return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]
        c_inner = [[R,0], [R,R*t], [R*t,R], [0,R]]
        c_outer = [[L,0], [L,L], [L,L], [0,L]]
        c_mid = [mid(c_inner[i], c_outer[i]) for i in range(4)]
        
        # Interleave CPs and weights properly
        cps_raw = numpy.array([c_inner, c_mid, c_outer]).transpose(1,0,2).reshape(-1, 2)
        weights_raw = numpy.array([[1,w_arc,w_arc,1], [1,1,1,1], [1,1,1,1]]).T.flatten()
        
        # Setup Symbolic Mapping
        W_func = function.dot(basis, weights_raw)
        X_func = function.dot(basis, cps_raw * weights_raw[:, None]) / W_func
        
        ns = Namespace()
        ns.x = X_func
        ns.N = basis
        ns.E = req.E
        ns.nu = req.nu
        ns.Tx = req.traction
        ns.R0 = R
        ns.L = L
        ns.fac = req.E / (1 - req.nu**2)
        ns.delta = numpy.eye(2)
        
        # 3. Physics
        ns.ubasis = basis.vector(2)
        ns.u = 'ubasis_ni ?lhs_n'
        ns.eps = '(u_i,j + u_j,i) / 2'
        ns.sigma = 'fac * ((1 - nu) * eps_ij + nu * eps_kk * delta_ij)'
        
        # 4. DBC & Traction
        # Symmetric boundary conditions
        sqr = topo.boundary['bottom'].integral('u_1^2' @ ns, degree=4)
        sqr += topo.boundary['left'].integral('u_0^2' @ ns, degree=4)
        cons = solver.optimize('lhs', sqr, droptol=1e-12)
        
        # Analytical Traction (Kirsch)
        ns.r = 'sqrt(x_i x_i)'
        ns.cost, ns.sint = 'x_0 / r', 'x_1 / r'
        ns.cos2t, ns.sin2t = 'cost^2 - sint^2', '2 cost sint'
        ns.cos4t, ns.sin4t = 'cos2t^2 - sin2t^2', '2 cos2t sin2t'
        ns.sxx_ex = 'Tx * (1 - (R0/r)^2 * (1.5 * cos2t + cos4t) + 1.5 * (R0/r)^4 * cos4t)'
        ns.syy_ex = '-Tx * ((R0/r)^2 * (0.5 * cos2t - cos4t) + 1.5 * (R0/r)^4 * cos4t)'
        ns.sxy_ex = '-Tx * ((R0/r)^2 * (0.5 * sin2t + sin4t) - 1.5 * (R0/r)^4 * sin4t)'
        ns.sigma_ex = function.stack([function.stack([ns.sxx_ex, ns.sxy_ex]), function.stack([ns.sxy_ex, ns.syy_ex])])
        ns.t = 'sigma_ex_ij n_j'
        
        res = topo.integral('sigma_ij ubasis_ni,j' @ ns, degree=4)
        res -= topo.boundary['right'].integral('t_i ubasis_ni' @ ns, degree=4)
        res -= topo.boundary['top'].integral('t_i ubasis_ni' @ ns, degree=4)
        
        # System Solution
        lhs = solver.solve_linear('lhs', res, constrain=cons)
        
        # 5. Result Sampling
        # Sample the stress field on a uniform grid for the web heatmap
        # res=30 matches JS heatmap resolution
        grid = topo.sample('bezier', 15) #bezier sampling for clean grid
        x_pts, sxx, syy, sxy = grid.eval(['x', 'sigma_00', 'sigma_11', 'sigma_01'] @ ns, lhs=lhs)
        
        return {
            "points": x_pts.tolist(), 
            "sxx": sxx.tolist(),
            "syy": syy.tolist(),
            "sxy": sxy.tolist(),
            "message": "Solved successfully with Nutils"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
