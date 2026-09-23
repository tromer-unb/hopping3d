"""Allow ``python -m hopping3d config.json`` as a portable CLI fallback."""

from .cli import main


if __name__ == "__main__":
    raise SystemExit(main())
