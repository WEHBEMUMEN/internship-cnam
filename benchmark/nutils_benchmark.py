import nutils, numpy
from nutils import mesh, function, solver, export, cli
from nutils.expression_v2 import Namespace

# Parameters matching Fig 4.2 / JS solver
R = 1.0; L = 4.0; E = 1e5; nu = 0.3; Tx = 100

def solve_kirsch(nrefine=0):
    # 1. Parameter Space
    topo, (u, v) = mesh.rectilinear([numpy.linspace(0, 1, 3), numpy.linspace(0, 1, 2)])
    for _ in range(nrefine): topo = topo.refine(1)
    
    # Quadratic Basis
    basis = topo.basis('spline', degree=2)
    
    # 2. Geometry (NURBS Fig 4.2)
    t = numpy.tan(numpy.pi/8); w_arc = numpy.cos(numpy.pi/8)
    def mid(p1, p2): return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]
    c_inner = [[R,0], [R,R*t], [R*t,R], [0,R]]
    c_outer = [[L,0], [L,L], [L,L], [0,L]]
    c_mid = [mid(c_inner[i], c_outer[i]) for i in range(4)]
    cps = numpy.array([c_inner, c_mid, c_outer]).transpose(1,0,2).reshape(-1, 2)
    weights = numpy.array([[1,w_arc,w_arc,1], [1,1,1,1], [1,1,1,1]]).T.flatten()
    
    # Use symbolic functions for mapping to avoid index collisions
    W_func = function.dot(basis, weights)
    X_func = function.dot(basis, cps * weights[:, None]) / W_func
    
    ns = Namespace()
    ns.x = X_func
    ns.N = basis
    ns.E = E; ns.nu = nu; ns.Tx = Tx; ns.R0 = R; ns.L = L
    ns.fac = E / (1 - nu**2)
    ns.delta = numpy.eye(2)
    
    # 3. Physics
    ns.ubasis = basis.vector(2)
    ns.u = 'ubasis_ni ?lhs_n'
    ns.eps = '(u_i,j + u_j,i) / 2'
    ns.sigma = 'fac * ((1 - nu) * eps_ij + nu * eps_kk * delta_ij)'
    
    # 4. Boundaries & Traction
    ns.r = 'sqrt(x_i x_i)'
    ns.cost, ns.sint = 'x_0 / r', 'x_1 / r'
    ns.cos2t, ns.sin2t = 'cost^2 - sint^2', '2 cost sint'
    ns.cos4t, ns.sin4t = 'cos2t^2 - sin2t^2', '2 cos2t sin2t'
    ns.sxx_ex = 'Tx * (1 - (R0/r)^2 * (1.5 * cos2t + cos4t) + 1.5 * (R0/r)^4 * cos4t)'
    ns.syy_ex = '-Tx * ((R0/r)^2 * (0.5 * cos2t - cos4t) + 1.5 * (R0/r)^4 * cos4t)'
    ns.sxy_ex = '-Tx * ((R0/r)^2 * (0.5 * sin2t + sin4t) - 1.5 * (R0/r)^4 * sin4t)'
    ns.sigma_ex = function.stack([function.stack([ns.sxx_ex, ns.sxy_ex]), function.stack([ns.sxy_ex, ns.syy_ex])])
    ns.t = 'sigma_ex_ij n_j'
    
    sqr = topo.boundary['bottom'].integral('u_1^2' @ ns, degree=4)
    sqr += topo.boundary['left'].integral('u_0^2' @ ns, degree=4)
    cons = solver.optimize('lhs', sqr, droptol=1e-12)
    
    res = topo.integral('sigma_ij ubasis_ni,j' @ ns, degree=4)
    res -= topo.boundary['right'].integral('t_i ubasis_ni' @ ns, degree=4)
    res -= topo.boundary['top'].integral('t_i ubasis_ni' @ ns, degree=4)
    
    lhs = solver.solve_linear('lhs', res, constrain=cons)
    
    # 5. Error
    ns.diff2 = '(sigma_ex_ij - sigma_ij) (sigma_ex_ij - sigma_ij)'
    ns.exact2 = 'sigma_ex_ij sigma_ex_ij'
    # Exclude singularity
    ns.r_max = 0.85 * L * numpy.sqrt(2)
    ns.mask = 'r < r_max'
    
    err = numpy.sqrt(topo.integral('diff2 mask' @ ns, degree=8).eval(lhs=lhs))
    val = numpy.sqrt(topo.integral('exact2 mask' @ ns, degree=8).eval(lhs=lhs))
    
    l2 = 100 * err / val
    print(f"Refine {nrefine}: L2 Stress Error = {l2:.2f}%")
    return l2

if __name__ == "__main__":
    for i in range(4): solve_kirsch(i)
