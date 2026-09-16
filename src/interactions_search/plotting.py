"""PNGs 3D de la envolvente convexa (convex hull) del sitio activo y de cada
pocket hidrofóbico calificado: dispersión de puntos + vértices, o superficie
sólida coloreada por altura."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')  # headless: sin esto matplotlib puede requerir un $DISPLAY inexistente en servidores
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

__all__ = ["plot_hull_volume", "plot_hull_surface", "plot_ramachandran", "plot_chi_profile"]


def plot_hull_volume(points, hull, title, filename):
    """PNG de dispersión 3D: todos los puntos (negro) y los vértices de su
    envolvente convexa (rojo), con el volumen en el título."""
    vertices = points[hull.vertices]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], marker='.', color='black')
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], marker='x', color='red')
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def plot_hull_surface(points, hull, title, filename):
    """PNG de superficie 3D sólida de la envolvente convexa: cada cara
    triangular del hull (hull.simplices) coloreada según su altura promedio
    (colormap viridis) — vista tipo malla sólida, complementaria a
    plot_hull_volume() (dispersión de puntos + vértices)."""
    faces = points[hull.simplices]  # (n_faces, 3, 3)
    face_z = faces[:, :, 2].mean(axis=1)
    norm = plt.Normalize(face_z.min(), face_z.max()) if face_z.max() > face_z.min() \
        else plt.Normalize(face_z.min() - 1, face_z.max() + 1)
    colors = plt.cm.viridis(norm(face_z))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    poly = Poly3DCollection(faces, facecolors=colors, edgecolor='k', linewidths=0.3, alpha=0.95)
    ax.add_collection3d(poly)

    # add_collection3d no autoescala los límites de los ejes: hay que fijarlos
    # a mano con el bounding box de los puntos, si no el plot queda vacío.
    mins, maxs = points.min(axis=0), points.max(axis=0)
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def plot_ramachandran(df_phi_psi, title, filename):
    """PNG de dispersión phi (X) vs psi (Y) para los residuos del sitio
    activo (df_phi_psi: columnas Pos, Residue, phi, psi — de
    ramachandran.compute_active_site_phi_psi). Filas con phi o psi en None
    (vecino faltante/gap de secuencia) se descartan. GLY (sin restricción de
    Cβ) y PRO (anillo que fija phi) se destacan aparte por ser los outliers
    esperados de un plot de Ramachandran estándar."""
    df = df_phi_psi.dropna(subset=['phi', 'psi'])
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.axhline(0, color='0.85', linewidth=0.8, zorder=0)
    ax.axvline(0, color='0.85', linewidth=0.8, zorder=0)
    ax.set_xlabel('phi (°)')
    ax.set_ylabel('psi (°)')
    ax.set_title(title)

    is_gly = df['Residue'] == 'GLY'
    is_pro = df['Residue'] == 'PRO'
    other  = ~(is_gly | is_pro)
    ax.scatter(df.loc[other, 'phi'], df.loc[other, 'psi'], marker='o', color='black', label='other')
    ax.scatter(df.loc[is_gly, 'phi'], df.loc[is_gly, 'psi'], marker='^', color='tab:green', label='GLY')
    ax.scatter(df.loc[is_pro, 'phi'], df.loc[is_pro, 'psi'], marker='s', color='tab:orange', label='PRO')
    for _, row in df.iterrows():
        ax.annotate(f"{row['Residue']}{row['Pos']}", (row['phi'], row['psi']),
                   fontsize=6, xytext=(3, 3), textcoords='offset points')
    if len(df):
        ax.legend(loc='upper right', fontsize=8)
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def plot_chi_profile(df_chi, chi_name, title, filename):
    """PNG de dispersión de un único ángulo chi (chi_name: 'chi1'..'chi5')
    para los residuos del sitio activo: eje X = residuos (uno por posición,
    en el orden en que aparecen en df_chi), eje Y = ese chi en grados
    (-180°, 180°], mismo rango y grillas que plot_ramachandran() para que
    ambos se lean con la misma escala. Filas sin ese chi (residuo sin ese
    ángulo, ej. ALA no tiene ninguno, VAL solo tiene chi1) se descartan."""
    df = df_chi.dropna(subset=[chi_name])
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(df) + 1), 5))
    ax.set_ylim(-180, 180)
    ax.set_yticks(range(-180, 181, 60))
    ax.axhline(0, color='0.85', linewidth=0.8, zorder=0)
    ax.set_ylabel(f'{chi_name} (°)')
    ax.set_title(title)

    labels = [f"{row['Residue']}{row['Pos']}" for _, row in df.iterrows()]
    ax.scatter(range(len(df)), df[chi_name], marker='o', color='black', zorder=3)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_xlim(-0.5, max(len(df) - 0.5, 0.5))
    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)
