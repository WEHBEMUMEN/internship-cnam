import nutils, numpy
from nutils import mesh, function, solver
from nutils.expression_v2 import Namespace

def solve(req):
    # 1. Reconstruct current topology from UI state
    nU = len(req.cps)
    nV = len(req.cps[0])
    
    # Define custom topology from knots
    topo, (u, v) = mesh.rectilinear([req.U, req.V])
    basis = topo.basis('spline', degree=(req.p, req.q))
    
    # 2. Geometry: Map Nutils local (u,v) to our NURBS physical (x,y)
    cps_flat = []
    ws_flat = []
    for i in range(nU):
        for j in range(nV):
            cps_flat.append([req.cps[i][j].x, req.cps[i][j].y])
            ws_flat.append(req.cps[i][j].w)
    
    cps_raw = numpy.array(cps_flat)
    weights_raw = numpy.array(ws_flat)
    
    # Mapping functions
    W_func = function.dot(basis, weights_raw)
    X_func = function.dot(basis, cps_raw * weights_raw[:, None]) / W_func
    
    ns = Namespace()
    ns.x = X_func
    ns.N = basis
    ns.E = req.E; ns.nu = req.nu; ns.Tx = req.traction
    ns.R0 = req.radius; ns.L = req.length
    ns.fac = req.E / (1 - req.nu**2)
    ns.delta = numpy.eye(2)
    
    # 3. Physics
    ns.ubasis = basis.vector(2)
    ns.u = 'ubasis_ni ?lhs_n'
    ns.eps = '(u_i,j + u_j,i) / 2'
    ns.sigma = 'fac * ((1 - nu) * eps_ij + nu * eps_kk * delta_ij)'
    
    # 4. DBCs: Standard Kirsch benchmark boundaries
    sqr = topo.boundary['bottom'].integral('u_1^2' @ ns, degree=4)
    sqr += topo.boundary['left'].integral('u_0^2' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-12)
    
    # 5. Analytical Traction
    ns.r = 'sqrt(x_i x_i)'
    ns.cost, ns.sint = 'x_0 / r', 'x_1 / r'
    ns.cos2t, ns.sint2t = 'cost^2 - sint^2', '2 cost sint' # Fix variable name from server.py
    ns.cos4t, ns.sin4t = 'cos2t^2 - sint2t^2', '2 cos2t sint2t'
    ns.sxx_ex = 'Tx * (1 - (R0/r)^2 * (1.5 * cos2t + cos4t) + 1.5 * (R0/r)^4 * cos4t)'
    ns.syy_ex = '-Tx * ((R0/r)^2 * (0.5 * cos2t - cos4t) + 1.5 * (R0/r)^4 * cos4t)'
    ns.sxy_ex = '-Tx * ((R0/r)^2 * (0.5 * sint2t + sin4t) - 1.5 * (R0/r)^4 * sin4t)'
    ns.sigma_ex = function.stack([function.stack([ns.sxx_ex, ns.sxy_ex]), function.stack([ns.sxy_ex, ns.syy_ex])])
    ns.t = 'sigma_ex_ij n_j'
    
    res = topo.integral('sigma_ij ubasis_ni,j' @ ns, degree=4)
    res -= topo.boundary['right'].integral('t_i ubasis_ni' @ ns, degree=4)
    res -= topo.boundary['top'].integral('t_i ubasis_ni' @ ns, degree=4)
    
    # 6. Solve
    lhs = solver.solve_linear('lhs', res, constrain=cons)
    return lhs.tolist()
