# AGENTS.md

## General principles

- Always write clear, maintainable, and deterministic code.
- Prefer simple and robust implementations over overly clever solutions.
- Keep responsibilities separated across scripts.
- Preserve stable input/output formats once they are defined in the spec.
- Do not silently change file formats, column names, argument names, or output ordering.
- When behavior is ambiguous, follow the spec files in the repository.
- If the spec and code disagree, update the code to match the spec unless explicitly instructed otherwise.

---

## Python rules

- Use Python 3.
- Write production-style scripts, not notebook-style code.
- Add functions for logically separate tasks.
- Include a `main()` entry point.
- Use `argparse` for CLI interfaces.
- Use explicit imports only.
- Avoid unnecessary dependencies.
- Prefer standard library, `numpy`, and `pandas` unless otherwise needed.
- Each Python script must include example commands showing how to run it, written as comments at the top of the file.
- Example commands should be written in triple-quoted docstrings (`""" ... """`) instead of `#` comments.
---

## Absolutely do not create pycache files

- Never create `__pycache__` directories.
- Never rely on `.pyc` files.
- Prevent bytecode generation in Python scripts.

Required practice:

- set `sys.dont_write_bytecode = True` in Python entry scripts
- avoid any workflow that generates `__pycache__`
- if needed, document running with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 script.py