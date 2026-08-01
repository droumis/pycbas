"""Command-line interface for CBAS.

Usage:
    cbas gui              # launch the interactive GUI (default port 5007)
    cbas gui --port 5008  # custom port
    cbas --help           # show help
"""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        prog="cbas",
        description="CBAS — Choice-Wide Behavioral Association Study",
    )
    subparsers = parser.add_subparsers(dest="command")

    gui_parser = subparsers.add_parser("gui", help="Launch the interactive GUI")
    gui_parser.add_argument("--port", type=int, default=5007, help="Port (default: 5007)")
    gui_parser.add_argument("--no-browser", action="store_true", help="Don't auto-open browser")

    args = parser.parse_args()

    if args.command == "gui":
        launch_gui(args)
    else:
        parser.print_help()


def launch_gui(args):
    try:
        import panel  # noqa: F401
    except ImportError:
        print("The GUI requires additional dependencies. Install with:")
        print("  pip install pycbas[gui]")
        print("  # or: pipx install pycbas[gui]")
        sys.exit(1)

    app_path = Path(__file__).parent / "_app.py"
    if not app_path.exists():
        print(f"Cannot find GUI app at: {app_path}")
        sys.exit(1)

    cmd = ["panel", "serve", str(app_path), "--autoreload", "--port", str(args.port)]
    if not args.no_browser:
        cmd.append("--show")

    print(f"Starting CBAS GUI on http://localhost:{args.port}")
    sys.exit(subprocess.call(cmd))
