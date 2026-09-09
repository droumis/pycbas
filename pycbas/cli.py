"""Command-line interface for pyCBAS.

Usage:
    pycbas gui               # launch the interactive GUI
    pycbas gui --port 5008   # custom port
    pycbas --version         # print the installed version
    pycbas --help            # show help

There is no analysis subcommand: an analysis is either run from Python, which the
guide covers, or from the GUI. `python -m pycbas` is equivalent to `pycbas`.
"""

import argparse
import socket
import subprocess
import sys
from pathlib import Path


def _find_open_port(start=5007, end=5099):
    for port in range(start, end + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("localhost", port)) != 0:
                return port
    return start


def main():
    from . import __version__

    parser = argparse.ArgumentParser(
        prog="pycbas",
        description="pyCBAS — Choice-Wide Behavioral Association Study",
    )
    # Worth having, because 0.2.0 changed numerical output at ties: "which version
    # produced this result" is a question users need to be able to answer.
    parser.add_argument("--version", action="version",
                        version=f"pycbas {__version__}")
    subparsers = parser.add_subparsers(dest="command")

    gui_parser = subparsers.add_parser("gui", help="Launch the interactive GUI")
    gui_parser.add_argument("--port", type=int, default=None, help="Port (auto-selects open port if not specified)")
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
        sys.exit(1)

    app_path = Path(__file__).parent / "_app.py"
    if not app_path.exists():
        print(f"Cannot find GUI app at: {app_path}")
        sys.exit(1)

    port = args.port if args.port else _find_open_port()

    # No --autoreload: it is a development flag, and panel emits a FutureWarning
    # about watchfiles that would otherwise be the first thing a new user sees.
    cmd = ["panel", "serve", str(app_path), "--port", str(port)]
    if not args.no_browser:
        cmd.append("--show")

    print(f"Starting pyCBAS GUI on http://localhost:{port}")
    sys.exit(subprocess.call(cmd))


if __name__ == "__main__":
    # So that `python -m pycbas.cli` works rather than doing nothing at all, which is
    # what it did before: the module defined main() and never called it.
    main()
