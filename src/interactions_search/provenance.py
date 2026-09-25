"""Per-analysis configuration snapshots and provenance, including failed runs."""
import hashlib
import inspect
import json
import platform
import subprocess
import time
import uuid
import warnings
from datetime import datetime, timezone
from functools import wraps
from importlib import metadata
from pathlib import Path

import yaml
from rdkit import Chem

from interactions_search.config import InteractionConfig
from interactions_search.interaction_rules import REASON_CODES, hbond_search_cutoff

_ALIASES = {
    'options': {
        'interaction_coord': 'Interaction_Coord_Source', 'volume_plot': 'Volume_Plot',
        'bias': 'Bias', 'bias_validated_only': 'Bias_Validated_Only',
    },
    'distancias': {'Distances_Hidrofobica': 'Distancia_Hidrofobica',
                   'centroid_distance': 'Distancia_Centro_Activo'},
    'pockets': {'min_residues': 'Pocket_Min_Residues',
                'coverage_threshold': 'Pocket_Coverage_Threshold',
                'density_radius': 'Pocket_Density_Radius'},
}


def configuration_snapshot(cfg, hp=None):
    """Convert the effective legacy API dict into a valid, reusable CLI YAML."""
    data = InteractionConfig(acceptors={}, donors={}, acceptors_antecedent={}).model_dump()
    for section in ('options', 'distancias', 'angulos', 'aromaticidad', 'pockets'):
        for field in data[section]:
            key = _ALIASES.get(section, {}).get(field, field)
            if key in cfg:
                data[section][field] = cfg[key]
    for field, key in (('acceptors', 'Aceptores_Prot'), ('donors', 'Dadores_Prot'),
                       ('acceptors_antecedent', 'Aceptot_antecedent'), ('special', 'Special_case')):
        data[field] = cfg.get(key, {})
    if hp is not None:
        data['hotspot_pocket'] = hp.model_dump()
    return InteractionConfig.model_validate(data).model_dump()


def _json_default(value):
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, 'tolist'):
        return value.tolist()
    raise TypeError(f'Cannot record provenance value of type {type(value).__name__}')


def _atomic_text(path, contents):
    temporary = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    try:
        temporary.write_text(contents, encoding='utf-8')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def file_record(path):
    path = Path(path).resolve()
    result = {'path': str(path), 'exists': path.is_file()}
    if result['exists']:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        result.update(sha256=digest.hexdigest(), size_bytes=path.stat().st_size)
    return result


def _software():
    package = Path(__file__).resolve().parent
    tree = hashlib.sha256()
    for path in sorted([*package.rglob('*.py'), *package.rglob('*.json')]):
        tree.update(path.relative_to(package).as_posix().encode() + b'\0')
        tree.update(path.read_bytes() + b'\0')
    versions = {}
    for name in ('interactions-search', 'rdkit', 'openbabel-wheel', 'biopython',
                 'numpy', 'pandas', 'scipy', 'pydantic', 'PyYAML', 'matplotlib'):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    result = {'python': platform.python_version(), 'dependencies': versions,
              'source_sha256': tree.hexdigest(), 'git_commit': None, 'git_dirty': None}
    try:
        root = package.parent.parent
        revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root,
                                  capture_output=True, text=True, timeout=3, check=True)
        status = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=all', '--',
                                 'src/interactions_search', 'pyproject.toml'], cwd=root,
                                capture_output=True, text=True, timeout=3, check=True)
        result.update(git_commit=revision.stdout.strip(), git_dirty=bool(status.stdout.strip()))
    except (OSError, subprocess.SubprocessError):
        pass  # An installed package need not be inside a Git checkout.
    return result


def _folder(mode, args):
    receptor = Path(args['receptor_pdb']).stem
    if mode == 'ligand':
        return Path(f"{receptor}_{Path(args['Ligand_imput']).stem}").resolve()
    if mode == 'probe':
        return Path(f"{receptor}_probe" + (f"_{args['label']}" if args['label'] else '')).resolve()
    if mode == 'site_bias':
        x, y, z = args['point']
        return Path(f'{receptor}_site_{x:.1f}_{y:.1f}_{z:.1f}').resolve()
    if mode == 'hotspots':
        return Path(f'{receptor}_hotspot_pockets').resolve()
    raise ValueError(f'Unknown analysis mode: {mode}')


class RunRecord:
    def __init__(self, mode, args):
        self.folder = _folder(mode, args)
        self.folder.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc)
        run_id = started.strftime('%Y%m%dT%H%M%S%fZ') + '_' + uuid.uuid4().hex[:8]
        self.history = self.folder / 'run_history' / run_id
        self.history.mkdir(parents=True)
        self.started = time.perf_counter()
        self.before = self._output_state()
        context = args.get('run_context') or {}
        self.sources = {'receptor': args['receptor_pdb']}
        if mode == 'ligand':
            self.sources['ligand'] = args['Ligand_imput']
        if mode == 'hotspots':
            root = Path(args['hotspot_dir'])
            for name in ('acceptors', 'donors', 'hydrophobics'):
                for path in sorted((root / name).glob('*')):
                    if path.name in ('clusters.csv', 'cluster_points.pdb') or path.suffix == '.dx':
                        self.sources[f'hotspots/{name}/{path.name}'] = path
        self.sources.update(context.get('source_files', {}))
        ignored = {'cfg', 'hp', 'ligand_reference', 'run_context'}
        self.data = {
            'schema_version': 1, 'interaction_csv_schema_version': 2,
            'run_id': run_id, 'mode': mode, 'status': 'running',
            'started_at_utc': started.isoformat(),
            'parameters': {k: v for k, v in args.items() if k not in ignored},
            'cli_arguments': context.get('cli_arguments'),
            'runtime_config': args['cfg'],
            'effective_hbond_search_distance_A': hbond_search_cutoff(args['cfg']),
            'inputs_before': {name: file_record(path) for name, path in self.sources.items()},
            'software': _software(), 'reason_codes': REASON_CODES,
        }
        snapshot = configuration_snapshot(args['cfg'], args.get('hp'))
        config_text = yaml.safe_dump(snapshot, sort_keys=False, allow_unicode=True)
        _atomic_text(self.folder / 'config_used.yml', config_text)
        _atomic_text(self.history / 'config_used.yml', config_text)
        self.data['config_snapshot'] = {
            'file': str((self.history / 'config_used.yml').relative_to(self.folder)),
            'sha256': file_record(self.history / 'config_used.yml')['sha256'],
        }
        reference = args.get('ligand_reference')
        if reference is not None:
            reference_path = self.history / 'ligand_reference.sdf'
            # Keep template atom order: a canonical SMILES alone could change
            # first-match tie breaking for chemically ambiguous symmetric groups.
            writer = Chem.SDWriter(str(reference_path))
            try:
                writer.write(Chem.Mol(reference))
            finally:
                writer.close()
            self.data['ligand_reference'] = {
                'file': str(reference_path.relative_to(self.folder)),
                'sha256': file_record(reference_path)['sha256'],
                'canonical_smiles': Chem.MolToSmiles(reference),
            }
        else:
            self.data['ligand_reference'] = None
        self._save()

    def _output_state(self):
        return {str(p.relative_to(self.folder)): (p.stat().st_size, p.stat().st_mtime_ns)
                for p in self.folder.rglob('*') if p.is_file()
                and 'run_history' not in p.relative_to(self.folder).parts
                and p.name not in ('run_metadata.json', 'config_used.yml')}

    def _save(self):
        contents = json.dumps(self.data, indent=2, ensure_ascii=False,
                              default=_json_default, allow_nan=False) + '\n'
        _atomic_text(self.history / 'run_metadata.json', contents)
        _atomic_text(self.folder / 'run_metadata.json', contents)

    def finish(self, error=None):
        status = 'completed' if error is None else 'failed'
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            status = 'interrupted'
        self.data.update(status=status, finished_at_utc=datetime.now(timezone.utc).isoformat(),
                         duration_seconds=round(time.perf_counter() - self.started, 6))
        if error is not None:
            self.data['error'] = {'type': type(error).__name__, 'message': str(error)}
        self.data['inputs_after'] = {name: file_record(path) for name, path in self.sources.items()}
        self.data['outputs_written'] = {
            name: file_record(self.folder / name)
            for name, stat in self._output_state().items() if self.before.get(name) != stat
        }
        self._save()


def record_analysis(mode):
    """Record all public analysis modes without changing their return values."""
    def decorate(function):
        signature = inspect.signature(function)

        @wraps(function)
        def run(*args, **kwargs):
            bound = signature.bind(*args, **kwargs)
            bound.apply_defaults()
            record = RunRecord(mode, bound.arguments)
            try:
                result = function(*args, **kwargs)
            except BaseException as error:
                try:
                    record.finish(error)
                except (OSError, TypeError, ValueError) as recording_error:
                    warnings.warn(f'Could not finish failed-run metadata: {recording_error}',
                                  RuntimeWarning, stacklevel=2)
                raise
            record.finish()
            return result

        return run
    return decorate
