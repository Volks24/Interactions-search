try:
    from Interactions_search import main, analyze_pair, carga_variables
except ModuleNotFoundError:
    import sys
    from pathlib import Path
    _repo_root = Path(__file__).parent.parent.parent
    if str(_repo_root) not in sys.path:
        sys.path.insert(0, str(_repo_root))
    from Interactions_search import main, analyze_pair, carga_variables

__all__ = ["main", "analyze_pair", "carga_variables"]
