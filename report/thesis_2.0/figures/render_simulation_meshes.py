import matplotlib.pyplot as plt
import numpy as np
import re

def parse_mesh_raw(filepath):
    with open(filepath, 'r') as f:
        text = f.read()
    
    # Extract coordinates and rgb colors
    # pattern: (x, y, z) [rgb={r,g,b}]
    pattern = re.compile(r'\(\s*([0-9.-]+)\s*,\s*([0-9.-]+)\s*,\s*([0-9.-]+)\s*\)\s*\[rgb=\{\s*([0-9.-]+)\s*,\s*([0-9.-]+)\s*,\s*([0-9.-]+)\s*\}\]')
    matches = pattern.findall(text)
    
    triangles = []
    colors = []
    
    for i in range(0, len(matches), 3):
        if i + 2 < len(matches):
            p0 = [float(matches[i][0]), float(matches[i][1])]
            p1 = [float(matches[i+1][0]), float(matches[i+1][1])]
            p2 = [float(matches[i+2][0]), float(matches[i+2][1])]
            c0 = [float(matches[i][3]), float(matches[i][4]), float(matches[i][5])]
            c1 = [float(matches[i+1][3]), float(matches[i+1][4]), float(matches[i+1][5])]
            c2 = [float(matches[i+2][3]), float(matches[i+2][4]), float(matches[i+2][5])]
            triangles.append([p0, p1, p2])
            colors.append([(c0[0]+c1[0]+c2[0])/3, (c0[1]+c1[1]+c2[1])/3, (c0[2]+c1[2]+c2[2])/3])
            
    # Mathematically scale the displacement
    # Assuming symmetry about y=0 for the undeformed beam,
    # the downward displacement U_y(x) is the mean Y at each X slice.
    pts = np.array([pt for tri in triangles for pt in tri])
    xs = np.round(pts[:,0], 4)
    unique_xs = np.unique(xs)
    
    # Compute U_y(x) for each unique X
    u_y = {xi: np.mean(pts[xs == xi, 1]) for xi in unique_xs}
    
    # Apply a scale factor to the displacement (e.g. 1.5x)
    scale_factor = 1.5
    for i in range(len(triangles)):
        for j in range(3):
            x_val = np.round(triangles[i][j][0], 4)
            uy = u_y[x_val]
            # Deformed Y = Undeformed Y + U_y
            # We want New Y = Undeformed Y + scale * U_y
            # New Y = (Y - U_y) + scale * U_y = Y + (scale - 1) * U_y
            triangles[i][j][1] += (scale_factor - 1.0) * uy

    return np.array(triangles), np.array(colors)

def render_figure(raw_path, out_pdf, is_stress=False):
    triangles, colors = parse_mesh_raw(raw_path)
    
    fig, ax = plt.subplots(figsize=(7.5, 2.2), dpi=300)
    
    from matplotlib.collections import PolyCollection
    poly = PolyCollection(triangles, facecolors=colors, edgecolors='black', linewidths=0.2, antialiased=True)
    ax.add_collection(poly)
    
    ax.set_xlim(-0.2, 6.7)
    ax.set_ylim(-1.25, 0.8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(out_pdf, bbox_inches='tight', pad_inches=0.01, transparent=True)
    plt.close()
    print(f"Rendered: {out_pdf}")

if __name__ == '__main__':
    render_figure('c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/chapters/5_simulation/stress_mesh_raw.tex',
                  'c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/figures/fom_stress_mesh.pdf', is_stress=True)
    render_figure('c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/chapters/5_simulation/displacement_mesh_raw.tex',
                  'c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/figures/fom_displacement_mesh.pdf', is_stress=False)
