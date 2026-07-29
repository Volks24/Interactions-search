"""PNGs 3D de la envolvente convexa (convex hull) del sitio activo y de cada
pocket hidrofóbico calificado: dispersión de puntos + vértices, o superficie
sólida coloreada por altura."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')  # headless: sin esto matplotlib puede requerir un $DISPLAY inexistente en servidores
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

__all__ = ["plot_hull_volume", "plot_hull_surface"]


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
