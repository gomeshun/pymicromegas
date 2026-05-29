# Repository Instructions

## Python Environment

- The Python environment is managed with `uv`.
- Run Python commands through `uv run <command>`, for example `uv run python script.py` or `uv run pytest`.
- Activate the virtual environment with `source .venv/bin/activate` only when an interactive shell inside the environment is required.
- Add dependencies with `uv add <package>`; do not run `pip install` directly.
- If `.venv` or `pyproject.toml` is missing, ask the user before creating or changing the environment.

## Language And File Creation

- Write repository-maintained instructions, documentation, and demo notebooks in English unless the user explicitly requests another language.
- Use English for Markdown cells and explanatory text in notebooks.
- Create new files with English filenames and English content by default.
- Keep code identifiers, comments, examples, and generated artifacts in English unless compatibility or domain conventions require otherwise.