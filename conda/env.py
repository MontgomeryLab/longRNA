#!/usr/bin/env python3
"""
Conda environment management for longRNA pipeline.

This script helps manage conda environment specifications and lock files.
Lock files ensure reproducible environments across different machines.

Usage:
    python env.py rebuild-locks --env-file environment.yml
    python env.py make-archive --env-name longRNA
"""

import platform
import subprocess
import argparse
import sys
import os

# Lock file names
osx_lockfile = "conda-osx-64.lock"
linux_lockfile = "conda-linux-64.lock"

# Development environment for conda-lock
dev_env_name = "longRNA_dev"
dev_deps = [
    "python=3.10",
    "conda-lock>=2.0",
    "conda-pack>=0.7.0"
]


def get_args():
    parser = argparse.ArgumentParser(
        description="Manage conda environments and lock files for longRNA"
    )
    subcommands = parser.add_subparsers(
        title="Available commands",
        required=True,
        dest="command"
    )

    # rebuild-locks command
    rebuild_locks = subcommands.add_parser(
        "rebuild-locks",
        help="Rebuild osx-64 and linux-64 lockfiles from environment file"
    )
    rebuild_locks.add_argument(
        "--env-file",
        default="environment.yml",
        help="The environment.yml file to use (default: environment.yml)"
    )

    # make-archive command
    make_archive = subcommands.add_parser(
        "make-archive",
        help="Create an environment archive for offline installation"
    )
    make_archive.add_argument(
        "--env-name",
        default="longRNA",
        help="The name of the environment to archive (default: longRNA)"
    )

    # install command
    install = subcommands.add_parser(
        "install",
        help="Install environment from lock file (for current platform)"
    )
    install.add_argument(
        "--name",
        default="longRNA",
        help="Name for the new environment (default: longRNA)"
    )

    return parser.parse_args()


def get_conda_envs():
    """Returns a dict of {env_name: env_path} for all conda environments."""
    proc = subprocess.run(
        ["conda", "env", "list", "--json"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8"
    )

    try:
        proc.check_returncode()
        import json
        envs = json.loads(proc.stdout)
        return {os.path.basename(p): p for p in envs.get("envs", [])}
    except subprocess.CalledProcessError:
        print(proc.stderr, file=sys.stderr)
        sys.exit("Failed to get conda environment list.")


def ensure_dev_env():
    """Ensure the development environment exists with conda-lock installed."""
    envs = get_conda_envs()

    if dev_env_name in envs:
        print(f"Dev environment '{dev_env_name}' found.")
        return

    print(f"Creating dev environment '{dev_env_name}'...")
    proc = subprocess.run(
        ["conda", "create",
         "--name", dev_env_name,
         "--channel", "conda-forge",
         "--yes"] + dev_deps,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8"
    )

    try:
        proc.check_returncode()
        print("Dev environment created.")
    except subprocess.CalledProcessError:
        print(proc.stdout)
        sys.exit("Failed to create dev environment.")


def rebuild_lockfiles(env_file):
    """Rebuilds the osx-64 and linux-64 lockfiles from the environment file."""

    ensure_dev_env()

    supported_platforms = ["osx-64", "linux-64"]
    env_ctx = ["conda", "run", "-n", dev_env_name]

    for plat in supported_platforms:
        lock_file = f"conda-{plat}.lock"
        lock_cmd = [
            "conda-lock", "lock",
            "--file", env_file,
            "--kind", "explicit",
            "--platform", plat,
            "--filename-template", "conda-{platform}.lock"
        ]

        print(f"Building {plat} lockfile...")
        proc = subprocess.run(
            env_ctx + lock_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8"
        )

        try:
            proc.check_returncode()
            print(f"  Created: {lock_file}")
        except subprocess.CalledProcessError:
            print(proc.stdout)
            sys.exit(f"Error building {plat} lockfile.")

    print("Lock files rebuilt successfully!")


def install_from_lock(env_name):
    """Install environment from the appropriate lock file for current platform."""

    # Detect platform
    system = platform.system()
    if system == "Darwin":
        lock_file = osx_lockfile
    elif system == "Linux":
        lock_file = linux_lockfile
    else:
        sys.exit(f"Unsupported platform: {system}")

    if not os.path.exists(lock_file):
        sys.exit(f"Lock file not found: {lock_file}\nRun 'python env.py rebuild-locks' first.")

    print(f"Installing from {lock_file}...")
    proc = subprocess.run(
        ["conda", "create",
         "--name", env_name,
         "--file", lock_file,
         "--yes"],
        encoding="utf-8"
    )

    if proc.returncode == 0:
        print(f"\nEnvironment '{env_name}' created successfully!")
        print(f"Activate with: conda activate {env_name}")
    else:
        sys.exit("Installation failed.")


def make_env_archive(env_name):
    """Create a portable archive of the environment using conda-pack."""

    ensure_dev_env()

    system = platform.system()
    os_name = {"Darwin": "osx", "Linux": "linux"}.get(system, system.lower())
    archive_name = f"{env_name}_{os_name}.tar.gz"

    env_ctx = ["conda", "run", "-n", dev_env_name]
    command = ["conda-pack", "-n", env_name, "-o", archive_name]

    print(f"Creating archive: {archive_name}...")
    proc = subprocess.run(env_ctx + command, encoding="utf-8")

    if proc.returncode == 0:
        print(f"Archive created: {archive_name}")
    else:
        sys.exit("Failed to create archive.")


def main():
    args = get_args()

    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    if args.command == "rebuild-locks":
        rebuild_lockfiles(args.env_file)
    elif args.command == "install":
        install_from_lock(args.name)
    elif args.command == "make-archive":
        make_env_archive(args.env_name)

    print("Done.")


if __name__ == "__main__":
    main()
