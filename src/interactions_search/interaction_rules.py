"""Shared candidate search and validation for real ligands and point probes.

Distances are strict (<); H-bond angles use min < angle <= max. Missing
geometry is rejected unless the caller explicitly requests distance-only
validation because the probe does not define an orientation.
"""
import numpy as np

from interactions_search.config import Angles, Distances
from interactions_search.geometry import angle_three_points

DEFAULT_DISTANCES = Distances()
DEFAULT_ANGLES = Angles()

REASON_CODES = {
    'distance_outside_cutoff': 'Distance does not satisfy the strict type-specific cutoff.',
    'angle_outside_range': 'Angle does not satisfy the configured range.',
    'missing_required_geometry': 'Required angle cannot be computed from available geometry.',
    'distance_and_angle_pass': 'Both distance and required angle satisfy the criteria.',
    'distance_pass': 'Distance-only criterion passed before output rounding/grouping.',
    'distance_only_no_probe_orientation':
        'Distance passed; a point probe cannot define this angle.',
    'steric_clash': 'Probe is closer to a receptor heavy atom than Probe_Clash_Distance.',
}

_DISTANCE_KEYS = {
    'acceptor': 'Distances_Hidrogen_Bonds',
    'donor': 'Distances_Hidrogen_Bonds',
    'aromatic': 'Distances_Aromatic',
    'hydrophobic': 'Distancia_Hidrofobica',
    'salt_bridge': 'Distances_Salt_Bridge',
    'pi_cation': 'Distances_Pi_Cation',
    'clash': 'Probe_Clash_Distance',
}


def distance_cutoff(kind, cfg):
    key = _DISTANCE_KEYS[kind]
    field = 'Distances_Hidrofobica' if kind == 'hydrophobic' else key
    return cfg.get(key, getattr(DEFAULT_DISTANCES, field))


def hbond_search_cutoff(cfg):
    return max(distance_cutoff('donor', cfg), cfg.get(
        'Hydrogen_Bond_Search_Distance', DEFAULT_DISTANCES.Hydrogen_Bond_Search_Distance))


def within_distance(distance, cutoff):
    """Scalar/array predicate, evaluated before output rounding."""
    values = np.asarray(distance)
    return np.isfinite(values) & (values >= 0) & (values < cutoff)


def neighbors(point, coordinates, cutoff):
    """Yield (positional index, raw distance) for candidates within cutoff."""
    coords = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    distances = np.linalg.norm(coords - np.asarray(point, dtype=float), axis=1)
    for index in np.flatnonzero(within_distance(distances, cutoff)):
        yield int(index), float(distances[index])


def atoms_near(point, site, atom_table, cutoff):
    """Shared receptor atom selection for hydrophobic and charged contacts."""
    mask = [atom in atom_table.get(residue, ())
            for residue, atom in zip(site['Residue'], site['Atom'])]
    selected = site.loc[np.asarray(mask, dtype=bool)]
    coordinates = selected[['X', 'Y', 'Z']].to_numpy(dtype=float)
    for index, distance in neighbors(point, coordinates, cutoff):
        yield selected.iloc[index], distance, coordinates[index]


def evaluate_hbond(distance, angle, cfg, *, angle_required=True):
    """Return (status, semicolon-separated reasons), retaining all failed criteria."""
    failures = []
    if not within_distance(distance, distance_cutoff('donor', cfg)):
        failures.append('distance_outside_cutoff')
    lower = cfg.get('Angle_Hidrogen_Bonds_Min', DEFAULT_ANGLES.Angle_Hidrogen_Bonds_Min)
    upper = cfg.get('Angle_Hidrogen_Bonds_Max', DEFAULT_ANGLES.Angle_Hidrogen_Bonds_Max)
    if not np.isfinite(angle):
        if angle_required:
            failures.append('missing_required_geometry')
    elif not lower < angle <= upper:
        failures.append('angle_outside_range')
    if failures:
        return 'No', ';'.join(failures)
    reason = ('distance_and_angle_pass' if np.isfinite(angle)
              else 'distance_only_no_probe_orientation')
    return 'Yes', reason


def hbond_status(distance, angle, cfg, *, angle_required=True):
    return evaluate_hbond(distance, angle, cfg, angle_required=angle_required)[0]


def evaluate_aromatic(distance, angle, cfg, *, angle_required=True):
    failures = []
    if not within_distance(distance, distance_cutoff('aromatic', cfg)):
        failures.append('distance_outside_cutoff')
    parallel = cfg.get('Aromatic_Parallel_Max', DEFAULT_ANGLES.Aromatic_Parallel_Max)
    tshaped = cfg.get('Aromatic_TShaped_Min', DEFAULT_ANGLES.Aromatic_TShaped_Min)
    if angle_required:
        if not np.isfinite(angle):
            failures.append('missing_required_geometry')
        elif not (0 <= angle <= 90 and (angle < parallel or angle > tshaped)):
            failures.append('angle_outside_range')
    if failures:
        return 'No', ';'.join(failures)
    return 'Yes', ('distance_and_angle_pass' if angle_required
                   else 'distance_only_no_probe_orientation')


def aromatic_status(distance, angle, cfg, *, angle_required=True):
    return evaluate_aromatic(distance, angle, cfg, angle_required=angle_required)[0]


def receptor_acceptor_angle(donor, acceptor, site, antecedents):
    """D-A-antecedent, common to ligand donors and donor probes.

    Missing or ambiguous receptor coordinates yield NaN (not distance-only).
    """
    name = antecedents.get(acceptor['Residue'], {}).get(acceptor['Atom'], 'C')
    antecedent = site[(site['Pos'] == acceptor['Pos']) & (site['Atom'] == name)]
    if len(antecedent) != 1:
        return np.nan
    return angle_three_points(donor, acceptor[['X', 'Y', 'Z']].to_numpy(dtype=float),
                              antecedent[['X', 'Y', 'Z']].to_numpy(dtype=float)[0])
