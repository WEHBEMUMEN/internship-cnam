import matplotlib.pyplot as plt
import numpy as np

# Set publication style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

# Colors matching the thesis theme (primaryblue, accentcyan, coral)
c_p1 = '#D9534F'  # Coral Red (Linear p=1)
c_p2 = '#1F5B9E'  # Primary Blue (Quadratic p=2)
c_p3 = '#008B8B'  # Teal/Cyan (Cubic p=3)

# Data
# DOFs for resolutions: 10x5, 20x10, 40x20, 80x40, 160x80
res_labels = ['10×5', '20×10', '40×20', '80×40', '160×80']

dofs_p1 = np.array([132, 462, 1722, 6642, 26082])
dofs_p2 = np.array([168, 528, 1848, 6888, 26568])
dofs_p3 = np.array([208, 598, 1978, 7138, 27058])

# Reference asymptotic displacement = 4.238 mm
u_ref = 4.238

# Tip displacements (mm) across dyadic resolutions
disp_p1 = np.array([3.120, 3.612, 3.895, 4.054, 4.142])
disp_p2 = np.array([3.538, 3.889, 4.083, 4.185, 4.238])
disp_p3 = np.array([3.742, 4.035, 4.192, 4.231, 4.238])

# Relative error (%) vs reference 4.238 mm
err_p1 = np.abs(u_ref - disp_p1) / u_ref * 100
err_p2 = np.abs(u_ref - disp_p2) / u_ref * 100
err_p3 = np.abs(u_ref - disp_p3) / u_ref * 100

# Left Plot: Tip Displacement Convergence
ax1.plot(dofs_p1, disp_p1, 'o-', color=c_p1, linewidth=2.2, markersize=7, label=r'Linear ($p=1, C^0$)')
ax1.plot(dofs_p2, disp_p2, 's-', color=c_p2, linewidth=2.2, markersize=7, label=r'Quadratic ($p=2, C^1$)')
ax1.plot(dofs_p3, disp_p3, '^-', color=c_p3, linewidth=2.2, markersize=7, label=r'Cubic ($p=3, C^2$)')
ax1.axhline(y=u_ref, color='gray', linestyle='--', linewidth=1.5, alpha=0.8, label=f'Reference ($U_y = {u_ref:.3f}$ mm)')

# Annotate baseline point on p=2
ax1.scatter([6888], [4.185], s=120, facecolors='none', edgecolors=c_p2, linewidths=2.5, zorder=5)
ax1.annotate(r'FOM Baseline ($80\times40$)', xy=(6888, 4.185), xytext=(2200, 3.75),
             arrowprops=dict(arrowstyle='->', color=c_p2, lw=1.5),
             fontsize=10, fontweight='bold', color=c_p2)

ax1.set_xscale('log')
ax1.set_xlabel('Degrees of Freedom (DOFs)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Vertical Tip Displacement $U_y$ (mm)', fontsize=12, fontweight='bold')
ax1.set_title('(a) Displacement Convergence under Static Load', fontsize=12, fontweight='bold', pad=10)
ax1.grid(True, which='both', linestyle=':', alpha=0.6)
ax1.legend(frameon=True, fontsize=10, loc='lower right')

# Right Plot: Relative Error log-log Convergence Rate
ax2.plot(dofs_p1, err_p1, 'o-', color=c_p1, linewidth=2.2, markersize=7, label=r'$k$-Refinement $p=1$ ($C^0$)')
ax2.plot(dofs_p2[:-1], err_p2[:-1], 's-', color=c_p2, linewidth=2.2, markersize=7, label=r'$k$-Refinement $p=2$ ($C^1$)')
ax2.plot(dofs_p3[:-1], err_p3[:-1], '^-', color=c_p3, linewidth=2.2, markersize=7, label=r'$k$-Refinement $p=3$ ($C^2$)')

# Annotate resolutions on p=2 curve
for dof, err, lbl in zip(dofs_p2[:-1], err_p2[:-1], res_labels[:-1]):
    offset_y = 1.25 if lbl != '80×40' else 0.7
    offset_x = 1.05 if lbl != '80×40' else 1.15
    ax2.annotate(f'{lbl}\n({err:.2f}%)', xy=(dof, err), xytext=(dof * offset_x, err * offset_y),
                 fontsize=8.5, color=c_p2, fontweight='semibold')

# Highlight baseline error
ax2.scatter([6888], [1.24], s=120, facecolors='none', edgecolors=c_p2, linewidths=2.5, zorder=5)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Degrees of Freedom (DOFs)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Relative Tip Displacement Error (%)', fontsize=12, fontweight='bold')
ax2.set_title(r'(b) Error Convergence across Refinement Orders ($k$-Refinement)', fontsize=12, fontweight='bold', pad=10)
ax2.grid(True, which='both', linestyle=':', alpha=0.6)
ax2.legend(frameon=True, fontsize=10, loc='upper right')

plt.tight_layout()
out_png = 'c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/figures/mesh_k_convergence.png'
out_pdf = 'c:/Users/Wehbe/Documents/internship-report/internship_report/thesis_2.0/figures/mesh_k_convergence.pdf'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
print(f'Successfully exported:\n  {out_png}\n  {out_pdf}')
