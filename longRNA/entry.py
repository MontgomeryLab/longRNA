#!/usr/bin/env python3
"""
longRNA - RNA-seq Analysis Pipeline CLI

Usage:
    longRNA salmonRun [--path PATH_FILE] [--config CONFIG_FILE]
    longRNA rsemRun [--path PATH_FILE] [--config CONFIG_FILE]

Examples:
    longRNA salmonRun                    # Run Salmon with default configs
    longRNA rsemRun                      # Run RSEM with default configs
    longRNA salmonRun --config custom.yml
    longRNA rsemRun --path custom_path.yml --config custom_config.yml
"""

import argparse
import subprocess
import sys
import time
import yaml
from pathlib import Path
from datetime import datetime


# Get script directory
script_dir = Path(__file__).parent.resolve()


def resolve_paths_recursive(data, base_dir):
    """Recursively resolve all relative paths in the config"""
    if isinstance(data, dict):
        # Handle CWL File object: {class: File, path: ...}
        if 'class' in data and data['class'] == 'File' and 'path' in data:
            resolved = data.copy()
            original_path = resolved['path']
            if not Path(original_path).is_absolute():
                abs_path = (base_dir / original_path).resolve()
                resolved['path'] = str(abs_path)
            return resolved
        else:
            # Recurse into dict
            return {k: resolve_paths_recursive(v, base_dir) for k, v in data.items()}
    elif isinstance(data, list):
        # Recurse into list
        return [resolve_paths_recursive(item, base_dir) for item in data]
    elif isinstance(data, str):
        # Plain string path (not a CWL File object)
        if not Path(data).is_absolute():
            try:
                # Only resolve if it looks like a path (contains / or ..)
                if '/' in data or data.startswith('..'):
                    return str((base_dir / data).resolve())
            except:
                pass
        return data
    else:
        return data


def merge_configs(path_config_file, run_config_file, merged_config_file, pipeline_name):
    """Merge path and config files into a processed config file"""

    # Check if config files exist
    if not path_config_file.exists():
        print(f"❌ Error: {path_config_file} not found!")
        print(f"   Please create the path config file.")
        sys.exit(1)

    if not run_config_file.exists():
        print(f"❌ Error: {run_config_file} not found!")
        print(f"   Please create the config file.")
        sys.exit(1)

    # Load both config files
    with open(path_config_file, 'r') as f:
        paths = yaml.safe_load(f)

    with open(run_config_file, 'r') as f:
        config = yaml.safe_load(f)

    # Resolve relative paths in the paths config (relative to path config file's directory)
    path_config_dir = path_config_file.parent if path_config_file.is_absolute() else Path.cwd()
    paths = resolve_paths_recursive(paths, path_config_dir)

    # Merge: config takes precedence for overlapping keys
    merged = {**paths, **config}

    # Add timestamp to run_directory to prevent overwriting previous results
    timestamp = datetime.now().strftime("%m-%d_%H_%M_%S")
    if 'run_directory' in merged:
        base_dir = merged['run_directory']
        parent = str(Path(base_dir).parent)
        base_name = Path(base_dir).name

        # Remove any existing timestamp suffix to get clean base name
        if base_name.startswith(f'{pipeline_name}_results'):
            base_clean = f'{pipeline_name}_results'
        else:
            base_clean = base_name

        if parent == '.':
            merged['run_directory'] = f"{base_clean}_{timestamp}"
        else:
            merged['run_directory'] = f"{parent}/{base_clean}_{timestamp}"
        print(f"📁 Output directory: {merged['run_directory']}")

    # Add timestamp for reference
    merged['run_date'] = datetime.now().strftime("%Y%m%d")
    merged['run_timestamp'] = timestamp

    # Save merged config
    with open(merged_config_file, 'w') as f:
        f.write(f"# Auto-generated merged configuration\n")
        f.write(f"# Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# DO NOT EDIT - Edit path and config files instead\n\n")
        yaml.dump(merged, f, default_flow_style=False, sort_keys=False)

    print(f"✅ Config files merged: {merged_config_file}")
    return merged, paths


def run_salmon(args):
    """Run the Salmon RNA-seq pipeline"""

    # Determine config files
    path_config_file = Path(args.path) if args.path else script_dir / "path_salmon.yml"
    run_config_file = Path(args.config) if args.config else script_dir / "config_salmon.yml"
    merged_config_file = script_dir / "processed_config_salmon.yml"

    print("=" * 70)
    print("🧬 Salmon RNA-seq Analysis Pipeline")
    print("=" * 70)
    print(f"Path config: {path_config_file}")
    print(f"Run config: {run_config_file}")
    print("Merging configuration files...")

    # Merge configs
    config_data, path_config_data = merge_configs(
        path_config_file, run_config_file, merged_config_file, 'salmon'
    )

    # Get workflow file
    workflow = script_dir / path_config_data.get('workflow_file', 'cwl/workflows/rnaseq_salmon_workflow.cwl')

    # Get temp directory
    tmp_dir_from_config = path_config_data.get('tmp_directory') or config_data.get('tmp_directory', '/home/shliu/cwl_tmp')
    if Path(tmp_dir_from_config).is_absolute():
        external_tmp = Path(tmp_dir_from_config)
    else:
        external_tmp = (script_dir / tmp_dir_from_config).resolve()
    external_tmp.mkdir(parents=True, exist_ok=True)

    # Get cache directory
    cache_dir_from_config = path_config_data.get('cache_directory')
    if cache_dir_from_config:
        if Path(cache_dir_from_config).is_absolute():
            cache_dir = Path(cache_dir_from_config)
        else:
            cache_dir = (script_dir / cache_dir_from_config).resolve()
    else:
        cache_dir = external_tmp.parent / "cwl_cache_salmon"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Command
    command = [
        "cwltool",
        "--leave-tmpdir",
        "--cachedir", str(cache_dir),
        "--timestamps",
        "--relax-path-checks",
        "--on-error", "continue",
        "--outdir", "../salmon_results",
        "--tmpdir-prefix", str(external_tmp) + "/",
        str(workflow),
        str(merged_config_file)
    ]

    print()
    print("=" * 70)
    print("Steps: fastp → salmon → DESeq2 analysis")
    print(f"Using temp directory: {external_tmp}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Start timing
    start_time = time.time()

    # Run and capture output
    result = subprocess.run(command, capture_output=True, text=True)

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print stderr
    if result.stderr:
        print(result.stderr, end='')

    # Print useful stdout lines (skip JSON)
    if result.stdout:
        lines = result.stdout.split('\n')
        json_depth = 0

        for line in lines:
            stripped = line.strip()

            # Track JSON depth
            if stripped.startswith('{'):
                json_depth += stripped.count('{')
                json_depth -= stripped.count('}')
                if json_depth > 0:
                    continue

            if json_depth > 0:
                json_depth += stripped.count('{')
                json_depth -= stripped.count('}')
                continue

            # Print non-JSON lines
            if stripped and not stripped.startswith('}'):
                print(line)

    # Format elapsed time
    elapsed_minutes = elapsed_time / 60
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)

    # Print summary
    print("=" * 70)
    if result.returncode == 0:
        print("✅ Pipeline completed successfully!")
        output_dir_name = config_data.get('run_directory', 'salmon_results')
        print(f"📁 Results in: {output_dir_name}/")
        print("   ├── fastp_cleaned/          (cleaned FASTQ files)")
        print("   ├── salmon_quants/          (gene quantifications)")
        print("   └── analysis_results/       (DESeq2 results & plots)")

        # Clean up temp directory on success
        print(f"\n🧹 Cleaning up temp files in {external_tmp}")
        import shutil
        try:
            shutil.rmtree(external_tmp)
            external_tmp.mkdir(parents=True, exist_ok=True)
            print("   Temp files cleaned!")
        except Exception as e:
            print(f"   Warning: Could not clean temp files: {e}")
    else:
        print(f"❌ Pipeline failed (exit code: {result.returncode})")
        print(f"💡 Tip: Temp files are in {external_tmp}")
        print(f"   You can delete them with: sudo rm -rf {external_tmp}/*")

    # Print timing information
    print(f"\n⏱️  Total runtime: {elapsed_minutes:.2f} minutes ({hours:02d}:{minutes:02d}:{seconds:02d})")
    print(f"   Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    sys.exit(result.returncode)


def run_rsem(args):
    """Run the RSEM + STAR RNA-seq pipeline"""

    # Determine config files
    path_config_file = Path(args.path) if args.path else script_dir / "path_rsem.yml"
    run_config_file = Path(args.config) if args.config else script_dir / "config_rsem.yml"
    merged_config_file = script_dir / "processed_config_rsem.yml"

    print("=" * 70)
    print("🧬 RSEM + STAR RNA-seq Analysis Pipeline")
    print("=" * 70)
    print(f"Path config: {path_config_file}")
    print(f"Run config: {run_config_file}")
    print("Merging configuration files...")

    # Merge configs
    config_data, path_config_data = merge_configs(
        path_config_file, run_config_file, merged_config_file, 'rsem'
    )

    # Get workflow file
    workflow = script_dir / path_config_data.get('workflow_file', 'cwl/workflows/rnaseq_rsem_workflow.cwl')

    # Get temp directory
    tmp_dir_from_config = path_config_data.get('tmp_directory') or config_data.get('tmp_directory', '../../cwl_tmp_rsem')
    if Path(tmp_dir_from_config).is_absolute():
        external_tmp = Path(tmp_dir_from_config)
    else:
        external_tmp = (script_dir / tmp_dir_from_config).resolve()
    external_tmp.mkdir(parents=True, exist_ok=True)

    # Get cache directory
    cache_dir_from_config = path_config_data.get('cache_directory')
    if cache_dir_from_config:
        if Path(cache_dir_from_config).is_absolute():
            cache_dir = Path(cache_dir_from_config)
        else:
            cache_dir = (script_dir / cache_dir_from_config).resolve()
    else:
        cache_dir = external_tmp.parent / "cwl_cache_rsem"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Command
    command = [
        "cwltool",
        "--leave-tmpdir",
        "--cachedir", str(cache_dir),
        "--timestamps",
        "--relax-path-checks",
        "--on-error", "continue",
        "--outdir", "../rsem_results",
        "--tmpdir-prefix", str(external_tmp) + "/",
        str(workflow),
        str(merged_config_file)
    ]

    print()
    print("=" * 70)
    print("Steps: fastp → RSEM (STAR) → DESeq2 analysis")
    print(f"Using temp directory: {external_tmp}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Start timing
    start_time = time.time()

    # Run and capture output
    result = subprocess.run(command, capture_output=True, text=True)

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print stderr
    if result.stderr:
        print(result.stderr, end='')

    # Print useful stdout lines (skip JSON)
    if result.stdout:
        lines = result.stdout.split('\n')
        json_depth = 0

        for line in lines:
            stripped = line.strip()

            # Track JSON depth
            if stripped.startswith('{'):
                json_depth += stripped.count('{')
                json_depth -= stripped.count('}')
                if json_depth > 0:
                    continue

            if json_depth > 0:
                json_depth += stripped.count('{')
                json_depth -= stripped.count('}')
                continue

            # Print non-JSON lines
            if stripped and not stripped.startswith('}'):
                print(line)

    # Format elapsed time
    elapsed_minutes = elapsed_time / 60
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)

    # Print summary
    print("=" * 70)
    if result.returncode == 0:
        print("✅ Pipeline completed successfully!")
        output_dir_name = config_data.get('run_directory', 'rsem_results')
        print(f"📁 Results in: {output_dir_name}/")
        print("   ├── fastp_cleaned/          (cleaned FASTQ files)")
        print("   ├── rsem_quants/            (gene quantifications)")
        print("   ├── rsem_bams/              (BAM alignment files)")
        print("   ├── rsem_logs/              (RSEM & STAR log files)")
        print("   ├── analysis_results_rsem/  (DESeq2 results & plots)")
        print("   └── gene_counts_rsem/       (raw & normalized counts)")

        # Clean up temp directory on success
        print(f"\n🧹 Cleaning up temp files in {external_tmp}")
        import shutil
        try:
            shutil.rmtree(external_tmp)
            external_tmp.mkdir(parents=True, exist_ok=True)
            print("   Temp files cleaned!")
        except Exception as e:
            print(f"   Warning: Could not clean temp files: {e}")
    else:
        print(f"❌ Pipeline failed (exit code: {result.returncode})")
        print(f"💡 Tip: Temp files are in {external_tmp}")
        print(f"   You can delete them with: sudo rm -rf {external_tmp}/*")

    # Print timing information
    print(f"\n⏱️  Total runtime: {elapsed_minutes:.2f} minutes ({hours:02d}:{minutes:02d}:{seconds:02d})")
    print(f"   Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    sys.exit(result.returncode)


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        prog='longRNA',
        description='RNA-seq Analysis Pipeline - Salmon and RSEM workflows',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s salmonRun                           # Run Salmon with default configs
  %(prog)s rsemRun                             # Run RSEM with default configs
  %(prog)s salmonRun --config custom.yml       # Use custom config file
  %(prog)s rsemRun --path path.yml --config config.yml
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Pipeline to run')

    # Salmon subcommand
    salmon_parser = subparsers.add_parser(
        'salmonRun',
        help='Run Salmon RNA-seq pipeline'
    )
    salmon_parser.add_argument(
        '--path',
        type=str,
        help='Path configuration file (default: path_salmon.yml)'
    )
    salmon_parser.add_argument(
        '--config',
        type=str,
        help='Pipeline configuration file (default: config_salmon.yml)'
    )
    salmon_parser.set_defaults(func=run_salmon)

    # RSEM subcommand
    rsem_parser = subparsers.add_parser(
        'rsemRun',
        help='Run RSEM + STAR RNA-seq pipeline'
    )
    rsem_parser.add_argument(
        '--path',
        type=str,
        help='Path configuration file (default: path_rsem.yml)'
    )
    rsem_parser.add_argument(
        '--config',
        type=str,
        help='Pipeline configuration file (default: config_rsem.yml)'
    )
    rsem_parser.set_defaults(func=run_rsem)

    # Parse arguments
    args = parser.parse_args()

    # If no command specified, print help
    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Run the appropriate function
    args.func(args)


if __name__ == '__main__':
    main()
