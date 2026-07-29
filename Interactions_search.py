"""Compatibility shim: kept at the repo root so `python Interactions_search.py ...`
still works exactly as before. All logic now lives in the `interactions_search`
package (src/interactions_search/); this file just re-exports the CLI entry point."""
from interactions_search.cli import main
from interactions_search.pipeline import analyze_pair, carga_variables

__all__ = ["main", "analyze_pair", "carga_variables"]

if __name__ == '__main__':
    main()
