#!/usr/bin/env python3
"""
long - RNA-seq Analysis Pipeline CLI

Usage:
    long [rsem] [--path PATH_FILE] [--config CONFIG_FILE]
    long salmon [--path PATH_FILE] [--config CONFIG_FILE]

Examples:
    long                                 # Run RSEM pipeline (default)
    long rsem                            # Run RSEM pipeline
    long salmon                          # Run Salmon pipeline
    long rsem --config custom.yml        # Use custom config file
    long salmon --path path.yml --config config.yml
"""

import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import time
import yaml
from pathlib import Path
from datetime import datetime


# Get script directory
script_dir = Path(__file__).parent.resolve()


def read_samples_csv(csv_path, fastq_directory):
    """
    Read samples CSV file and extract sample information.

    CSV format:
        Forward Reads,Reverse Reads,Sample/Group Name,Replicate Number,Control
        mRNA324_1.fq.gz,mRNA324_2.fq.gz,N2,1,TRUE

    Returns:
        dict with keys: in_fastq_r1, in_fastq_r2, sample_basenames, sample_groups, sample_csv_file
    """
    if not os.path.exists(csv_path):
        print(f"❌ Error: Samples CSV file not found: {csv_path}")
        sys.exit(1)

    in_fastq_r1 = []
    in_fastq_r2 = []
    sample_basenames = []
    sample_groups = []

    print(f"📋 Reading samples from: {csv_path}")

    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)

        for row in reader:
            # Get file paths
            r1_file = row.get('Forward Reads', '').strip()
            r2_file = row.get('Reverse Reads', '').strip()
            group_name = row.get('Sample/Group Name', '').strip()

            if not r1_file:
                continue  # Skip empty rows

            # Build full paths
            r1_path = os.path.join(fastq_directory, r1_file)
            in_fastq_r1.append({'class': 'File', 'path': r1_path})

            if r2_file:  # Paired-end
                r2_path = os.path.join(fastq_directory, r2_file)
                in_fastq_r2.append({'class': 'File', 'path': r2_path})

            # Extract sample basename (remove _1.fq.gz, _R1.fastq.gz, etc.)
            basename = r1_file
            for suffix in ['_1.fq.gz', '_1.fastq.gz', '_R1.fq.gz', '_R1.fastq.gz',
                          '_1.fq', '_1.fastq', '_R1.fq', '_R1.fastq']:
                if basename.endswith(suffix):
                    basename = basename[:-len(suffix)]
                    break
            sample_basenames.append(basename)
            sample_groups.append(group_name)

    print(f"   Found {len(sample_basenames)} samples: {', '.join(sample_basenames)}")
    print(f"   Groups: {', '.join(sorted(set(sample_groups)))}")

    result = {
        'in_fastq_r1': in_fastq_r1,
        'sample_basenames': sample_basenames,
        'sample_groups': sample_groups,
        'sample_csv_file': {'class': 'File', 'path': csv_path}
    }

    # Only add R2 if we have paired-end data
    if in_fastq_r2:
        result['in_fastq_r2'] = in_fastq_r2

    return result


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

    # If samples_csv is specified, read sample info from CSV file
    if 'samples_csv' in paths and paths['samples_csv']:
        samples_csv_path = paths['samples_csv']
        fastq_directory = paths.get('fastq_directory', os.path.dirname(samples_csv_path))

        # Read samples from CSV and add to paths
        sample_data = read_samples_csv(samples_csv_path, fastq_directory)
        paths.update(sample_data)

    # Merge: config takes precedence for overlapping keys
    merged = {**paths, **config}

    # Add timestamp to run_directory to prevent overwriting previous results
    timestamp = datetime.now().strftime("%y_%m-%d_%H_%M")
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
    """Run the Salmon RNA-seq pipeline with sample-by-sample processing"""
    # Determine config files
    path_config_file = Path(args.path) if args.path else script_dir / "path_salmon.yml"
    run_config_file = Path(args.config) if args.config else script_dir / "config_salmon.yml"
    merged_config_file = script_dir / "processed_config_salmon.yml"

    print("=" * 70)
    print("🧬 Salmon RNA-seq Analysis Pipeline (Sequential Mode)")
    print("=" * 70)
    print(f"Path config: {path_config_file}")
    print(f"Run config: {run_config_file}")
    print("Merging configuration files...")

    # Merge configs
    config_data, path_config_data = merge_configs(
        path_config_file, run_config_file, merged_config_file, 'salmon'
    )

    # Get directories
    tmp_dir_from_config = path_config_data.get('tmp_directory') or config_data.get('tmp_directory', '/home/shliu/cwl_tmp')
    if Path(tmp_dir_from_config).is_absolute():
        external_tmp = Path(tmp_dir_from_config)
    else:
        external_tmp = (script_dir / tmp_dir_from_config).resolve()
    external_tmp.mkdir(parents=True, exist_ok=True)

    cache_dir_from_config = path_config_data.get('cache_directory')
    if cache_dir_from_config:
        if Path(cache_dir_from_config).is_absolute():
            cache_dir = Path(cache_dir_from_config)
        else:
            cache_dir = (script_dir / cache_dir_from_config).resolve()
    else:
        cache_dir = external_tmp.parent / "cwl_cache_salmon"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Get output directory
    run_directory = Path(config_data.get('run_directory', '../salmon_results'))
    if not run_directory.is_absolute():
        run_directory = (script_dir / run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)

    # Create output subdirectories
    dir_fastp = run_directory / config_data.get('dir_name_fastp', 'fastp_cleaned')
    dir_salmon = run_directory / config_data.get('dir_name_salmon', 'salmon_quants')
    dir_analysis = run_directory / config_data.get('dir_name_analysis', 'analysis_results_salmon')
    dir_counts = run_directory / config_data.get('dir_name_gene_counts', 'gene_counts_salmon')

    for d in [dir_fastp, dir_salmon, dir_analysis, dir_counts]:
        d.mkdir(parents=True, exist_ok=True)

    # Get sample information
    sample_names = config_data.get('sample_basenames', [])
    sample_groups = config_data.get('sample_groups', [])
    fastq_r1_list = config_data.get('in_fastq_r1', [])
    fastq_r2_list = config_data.get('in_fastq_r2', [])
    is_paired = bool(fastq_r2_list)

    print()
    print("=" * 70)
    print(f"Processing {len(sample_names)} samples sequentially (with cleanup)")
    print(f"Output directory: {run_directory}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    start_time = time.time()

    # ========== STEP 1: Build Salmon Index (once) ==========
    print("\n" + "=" * 70)
    print("Step 1: Building Salmon Index")
    print("=" * 70)

    index_workflow = script_dir / "cwl/tools/salmon-index.cwl"
    index_config = {
        'transcriptome': config_data.get('transcriptome'),
        'threads': config_data.get('salmon_threads', config_data.get('threads', 8))
    }

    index_config_file = external_tmp / "index_config.yml"
    with open(index_config_file, 'w') as f:
        yaml.dump(index_config, f)

    index_outdir = external_tmp / "salmon_index_output"
    index_outdir.mkdir(parents=True, exist_ok=True)

    index_command = [
        "cwltool",
        "--cachedir", str(cache_dir),
        "--outdir", str(index_outdir),
        "--tmpdir-prefix", str(external_tmp) + "/",
        str(index_workflow),
        str(index_config_file)
    ]

    print(f"Building Salmon index...")
    index_result = subprocess.run(index_command, capture_output=True, text=True)

    if index_result.returncode != 0:
        print(f"❌ Index build failed!")
        print(index_result.stderr)
        sys.exit(1)

    # Find the index directory in output
    index_dir = None
    for item in index_outdir.iterdir():
        if item.is_dir():
            index_dir = item
            break

    if not index_dir:
        # Try to find in cache
        for cache_entry in cache_dir.iterdir():
            if cache_entry.is_dir():
                for f in cache_entry.iterdir():
                    if f.is_dir() and (f / "seq.bin").exists():
                        index_dir = f
                        break
                if index_dir:
                    break

    if not index_dir:
        print("❌ Could not find Salmon index directory!")
        sys.exit(1)

    print(f"✅ Index built: {index_dir}")

    # ========== STEP 2: Process Each Sample ==========
    print("\n" + "=" * 70)
    print("Step 2: Processing Samples (fastp → Salmon)")
    print("=" * 70)

    per_sample_workflow = script_dir / "cwl/workflows/salmon_per_sample.cwl"
    quant_dirs = []
    failed_samples = []

    for i, sample_name in enumerate(sample_names):
        print(f"\n--- Sample {i+1}/{len(sample_names)}: {sample_name} ---")
        sample_start = time.time()

        # Prepare sample config
        sample_config = {
            'fastq_r1': fastq_r1_list[i],
            'sample_name': sample_name,
            'index_dir': {'class': 'Directory', 'path': str(index_dir)},
            'library_type': config_data.get('library_type', 'A'),
            'threads': config_data.get('salmon_threads', config_data.get('threads', 8)),
        }

        if is_paired and i < len(fastq_r2_list):
            sample_config['fastq_r2'] = fastq_r2_list[i]

        # Add optional parameters
        for param in ['qualified_quality_phred', 'length_required', 'detect_adapter_for_pe',
                      'trim_poly_x', 'trim_poly_g']:
            if param in config_data and config_data[param] is not None:
                sample_config[param] = config_data[param]

        sample_config_file = external_tmp / f"sample_{sample_name}_config.yml"
        with open(sample_config_file, 'w') as f:
            yaml.dump(sample_config, f)

        sample_outdir = external_tmp / f"sample_{sample_name}_output"
        sample_outdir.mkdir(parents=True, exist_ok=True)

        # Run per-sample workflow
        sample_command = [
            "cwltool",
            "--cachedir", str(cache_dir),
            "--outdir", str(sample_outdir),
            "--tmpdir-prefix", str(external_tmp) + "/",
            str(per_sample_workflow),
            str(sample_config_file)
        ]

        print(f"   Running fastp + Salmon...")
        sample_result = subprocess.run(sample_command, capture_output=True, text=True)

        if sample_result.returncode != 0:
            print(f"   ❌ Sample {sample_name} failed!")
            if sample_result.stderr:
                err_lines = sample_result.stderr.strip().split('\n')
                for line in err_lines[-5:]:
                    print(f"      {line}")
            failed_samples.append(sample_name)
            continue

        # Move outputs to final directories
        print(f"   Moving outputs to final directory...")

        # Find and move files
        for f in sample_outdir.rglob('*'):
            if not f.is_file():
                continue

            fname = f.name

            # Fastp outputs
            if fname.endswith('_R1_clean.fastq.gz') or fname.endswith('_R2_clean.fastq.gz'):
                copy_file_to_dir(f, dir_fastp)
            elif fname.endswith('.html') and 'fastp' in str(f):
                copy_file_to_dir(f, dir_fastp)
            elif fname.endswith('.json') and 'fastp' in str(f):
                copy_file_to_dir(f, dir_fastp)

        # Copy quant directory
        for d in sample_outdir.rglob('*'):
            if d.is_dir() and (d / 'quant.sf').exists():
                dest_quant_dir = dir_salmon / sample_name
                if dest_quant_dir.exists():
                    shutil.rmtree(dest_quant_dir)
                shutil.copytree(d, dest_quant_dir)
                quant_dirs.append(str(dest_quant_dir))
                print(f"      {sample_name}/ → {dir_salmon.name}/")
                break

        # Clean up cache for this sample
        print(f"   Cleaning up cache...")
        deleted = cleanup_sample_cache(cache_dir, sample_name)
        print(f"      Deleted {deleted} cache entries")

        # Clean up sample output directory
        shutil.rmtree(sample_outdir, ignore_errors=True)

        sample_elapsed = time.time() - sample_start
        print(f"   ✅ {sample_name} complete ({sample_elapsed/60:.1f} min)")

    # ========== STEP 3: Run DESeq2 ==========
    if len(quant_dirs) > 0:
        print("\n" + "=" * 70)
        print("Step 3: Running DESeq2 Analysis")
        print("=" * 70)

        # Filter to only successfully processed samples
        successful_samples = [s for s in sample_names if s not in failed_samples]
        successful_groups = [sample_groups[i] for i, s in enumerate(sample_names) if s not in failed_samples]

        deseq_workflow = script_dir / "cwl/tools/deseq-analysis-salmon-standalone.cwl"
        deseq_config = {
            'r_script': {'class': 'File', 'path': str(script_dir / "scripts/deseq_analysis_salmon.R")},
            'gtf_annotation': config_data.get('gtf_annotation'),
            'quant_dirs': [{'class': 'Directory', 'path': p} for p in quant_dirs],
            'sample_names': successful_samples,
            'sample_groups': successful_groups,
            'output_prefix': config_data.get('run_name', 'deseq2_salmon')
        }

        # Add optional inputs
        if config_data.get('sample_csv_file'):
            deseq_config['sample_csv_file'] = config_data.get('sample_csv_file')
        if config_data.get('geneid2name_file'):
            deseq_config['geneid2name_file'] = config_data.get('geneid2name_file')

        deseq_config_file = external_tmp / "deseq_config.yml"
        with open(deseq_config_file, 'w') as f:
            yaml.dump(deseq_config, f)

        deseq_outdir = external_tmp / "deseq_output"
        deseq_outdir.mkdir(parents=True, exist_ok=True)

        deseq_command = [
            "cwltool",
            "--outdir", str(deseq_outdir),
            "--tmpdir-prefix", str(external_tmp) + "/",
            str(deseq_workflow),
            str(deseq_config_file)
        ]

        print(f"Running DESeq2 on {len(successful_samples)} samples...")
        deseq_result = subprocess.run(deseq_command, capture_output=True, text=True)

        if deseq_result.returncode == 0:
            print("✅ DESeq2 analysis complete!")

            # Move DESeq2 outputs
            print("   Moving DESeq2 outputs...")
            for f in deseq_outdir.rglob('*'):
                if not f.is_file():
                    continue
                fname = f.name
                if fname.endswith('.tsv'):
                    copy_file_to_dir(f, dir_counts, verbose=True)
                elif fname.endswith('.csv') or fname.endswith('.pdf'):
                    copy_file_to_dir(f, dir_analysis, verbose=True)
        else:
            print("❌ DESeq2 analysis failed!")
            if deseq_result.stderr:
                err_lines = deseq_result.stderr.strip().split('\n')
                for line in err_lines[-10:]:
                    print(f"   {line}")

    # ========== Final Summary ==========
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_minutes = elapsed_time / 60
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)

    print("\n" + "=" * 70)
    if len(failed_samples) == 0:
        print("✅ Pipeline completed successfully!")
    else:
        print(f"⚠️  Pipeline completed with {len(failed_samples)} failed samples:")
        for s in failed_samples:
            print(f"   - {s}")

    print(f"\n📁 Results in: {run_directory}/")
    print(f"   ├── {dir_fastp.name}/          (cleaned FASTQ files)")
    print(f"   ├── {dir_salmon.name}/         (Salmon quantifications)")
    print(f"   ├── {dir_analysis.name}/  (DESeq2 results & plots)")
    print(f"   └── {dir_counts.name}/      (raw & normalized counts)")

    print(f"\n⏱️  Total runtime: {elapsed_minutes:.2f} minutes ({hours:02d}:{minutes:02d}:{seconds:02d})")
    print(f"   Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    sys.exit(0 if len(failed_samples) == 0 else 1)


def move_file_to_dir(src_path, dest_dir, verbose=True):
    """Move a file to destination directory, creating dir if needed"""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    src = Path(src_path)
    if src.exists():
        dest = dest_dir / src.name
        shutil.move(str(src), str(dest))
        if verbose:
            print(f"      {src.name} → {dest_dir.name}/")
        return dest
    return None


def copy_file_to_dir(src_path, dest_dir, verbose=True):
    """Copy a file to destination directory, creating dir if needed"""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    src = Path(src_path)
    if src.exists():
        dest = dest_dir / src.name
        shutil.copy2(str(src), str(dest))
        if verbose:
            print(f"      {src.name} → {dest_dir.name}/")
        return dest
    return None


def find_files_in_cache(cache_dir, pattern):
    """Find files matching pattern in cache directory"""
    import fnmatch
    matches = []
    for root, dirs, files in os.walk(cache_dir):
        for f in files:
            if fnmatch.fnmatch(f, pattern):
                matches.append(Path(root) / f)
    return matches


def cleanup_sample_cache(cache_dir, sample_name, keep_genes_results=False):
    """Delete cache entries for a specific sample"""
    cache_path = Path(cache_dir)
    deleted_count = 0

    for cache_entry in cache_path.iterdir():
        if not cache_entry.is_dir():
            continue

        # Check if this cache entry contains files for this sample
        has_sample_files = False
        has_genes_results = False

        for f in cache_entry.rglob('*'):
            if sample_name in f.name:
                has_sample_files = True
                if f.name.endswith('.genes.results'):
                    has_genes_results = True

        # Delete if it has sample files (but keep if it only has genes.results and we want to keep those)
        if has_sample_files:
            if keep_genes_results and has_genes_results:
                continue
            try:
                shutil.rmtree(cache_entry)
                deleted_count += 1
            except Exception as e:
                print(f"      Warning: Could not delete {cache_entry}: {e}")

    return deleted_count


def run_rsem(args):
    """Run the RSEM + STAR RNA-seq pipeline with sample-by-sample processing"""
    # Determine config files
    path_config_file = Path(args.path) if args.path else script_dir / "path_rsem.yml"
    run_config_file = Path(args.config) if args.config else script_dir / "config_rsem.yml"
    merged_config_file = script_dir / "processed_config_rsem.yml"

    print("=" * 70)
    print("🧬 RSEM + STAR RNA-seq Analysis Pipeline (Sequential Mode)")
    print("=" * 70)
    print(f"Path config: {path_config_file}")
    print(f"Run config: {run_config_file}")
    print("Merging configuration files...")

    # Merge configs
    config_data, path_config_data = merge_configs(
        path_config_file, run_config_file, merged_config_file, 'rsem'
    )

    # Get directories
    tmp_dir_from_config = path_config_data.get('tmp_directory') or config_data.get('tmp_directory', '../../cwl_tmp_rsem')
    if Path(tmp_dir_from_config).is_absolute():
        external_tmp = Path(tmp_dir_from_config)
    else:
        external_tmp = (script_dir / tmp_dir_from_config).resolve()
    external_tmp.mkdir(parents=True, exist_ok=True)

    cache_dir_from_config = path_config_data.get('cache_directory')
    if cache_dir_from_config:
        if Path(cache_dir_from_config).is_absolute():
            cache_dir = Path(cache_dir_from_config)
        else:
            cache_dir = (script_dir / cache_dir_from_config).resolve()
    else:
        cache_dir = external_tmp.parent / "cwl_cache_rsem"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Get output directory
    run_directory = Path(config_data.get('run_directory', '../rsem_results'))
    if not run_directory.is_absolute():
        run_directory = (script_dir / run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)

    # Create output subdirectories
    dir_fastp = run_directory / config_data.get('dir_name_fastp', 'fastp_cleaned')
    dir_rsem = run_directory / config_data.get('dir_name_rsem', 'rsem_quants')
    dir_bams = run_directory / config_data.get('dir_name_rsem_bams', 'sorted_bam')
    dir_logs = run_directory / config_data.get('dir_name_rsem_logs', 'rsem_logs')
    dir_analysis = run_directory / config_data.get('dir_name_analysis', 'analysis_results_rsem')
    dir_counts = run_directory / config_data.get('dir_name_gene_counts', 'gene_counts_rsem')

    for d in [dir_fastp, dir_rsem, dir_bams, dir_logs, dir_analysis, dir_counts]:
        d.mkdir(parents=True, exist_ok=True)

    # Get sample information
    sample_names = config_data.get('sample_basenames', [])
    sample_groups = config_data.get('sample_groups', [])
    fastq_r1_list = config_data.get('in_fastq_r1', [])
    fastq_r2_list = config_data.get('in_fastq_r2', [])
    is_paired = bool(fastq_r2_list)

    print()
    print("=" * 70)
    print(f"Processing {len(sample_names)} samples sequentially (with cleanup)")
    print(f"Output directory: {run_directory}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    start_time = time.time()

    # ========== STEP 1: Build RSEM Reference (once) ==========
    print("\n" + "=" * 70)
    print("Step 1: Building RSEM + STAR Reference")
    print("=" * 70)

    ref_workflow = script_dir / "cwl/tools/build-rsem-star-reference.cwl"
    ref_config = {
        'genome_fasta': config_data.get('genome_fasta'),
        'gtf_annotation': config_data.get('gtf_annotation'),
        'reference_name': config_data.get('rsem_reference_name', 'rsem_reference'),
        'sjdboverhang': config_data.get('rsem_sjdboverhang', 100),
        'threads': config_data.get('threads', 8)
    }
    if config_data.get('rsem_genomeSAindexNbases'):
        ref_config['genomeSAindexNbases'] = config_data.get('rsem_genomeSAindexNbases')

    ref_config_file = external_tmp / "ref_config.yml"
    with open(ref_config_file, 'w') as f:
        yaml.dump(ref_config, f)

    ref_outdir = external_tmp / "rsem_ref_output"
    ref_outdir.mkdir(parents=True, exist_ok=True)

    ref_command = [
        "cwltool",
        "--cachedir", str(cache_dir),
        "--outdir", str(ref_outdir),
        "--tmpdir-prefix", str(external_tmp) + "/",
        str(ref_workflow),
        str(ref_config_file)
    ]

    print(f"Building reference (this may take a while)...")
    ref_result = subprocess.run(ref_command, capture_output=True, text=True)

    if ref_result.returncode != 0:
        print(f"❌ Reference build failed!")
        print(ref_result.stderr)
        sys.exit(1)

    # Find the reference directory in output
    reference_dir = None
    for item in ref_outdir.iterdir():
        if item.is_dir():
            reference_dir = item
            break

    if not reference_dir:
        # Try to find in cache
        for cache_entry in cache_dir.iterdir():
            if cache_entry.is_dir():
                for f in cache_entry.iterdir():
                    if f.is_dir() and (f / "SAindex").exists():
                        reference_dir = f
                        break
                if reference_dir:
                    break

    if not reference_dir:
        print("❌ Could not find RSEM reference directory!")
        sys.exit(1)

    print(f"✅ Reference built: {reference_dir}")

    # ========== STEP 2: Process Each Sample ==========
    print("\n" + "=" * 70)
    print("Step 2: Processing Samples (fastp → RSEM)")
    print("=" * 70)

    per_sample_workflow = script_dir / "cwl/workflows/rsem_per_sample.cwl"
    genes_results_files = []
    failed_samples = []

    for i, sample_name in enumerate(sample_names):
        print(f"\n--- Sample {i+1}/{len(sample_names)}: {sample_name} ---")
        sample_start = time.time()

        # Prepare sample config
        sample_config = {
            'fastq_r1': fastq_r1_list[i],
            'sample_name': sample_name,
            'is_paired_end': is_paired,
            'reference_dir': {'class': 'Directory', 'path': str(reference_dir)},
            'reference_name': config_data.get('rsem_reference_name', 'rsem_reference'),
            'forward_prob': config_data.get('rsem_forward_prob', 0.5),
            'threads': config_data.get('rsem_threads', config_data.get('threads', 8)),
        }

        if is_paired and i < len(fastq_r2_list):
            sample_config['fastq_r2'] = fastq_r2_list[i]

        # Add optional parameters
        for param in ['qualified_quality_phred', 'length_required', 'detect_adapter_for_pe',
                      'trim_poly_x', 'trim_poly_g', 'perfect_match_only', 'unique_mappers_only']:
            if param in config_data and config_data[param] is not None:
                sample_config[param] = config_data[param]

        sample_config_file = external_tmp / f"sample_{sample_name}_config.yml"
        with open(sample_config_file, 'w') as f:
            yaml.dump(sample_config, f)

        sample_outdir = external_tmp / f"sample_{sample_name}_output"
        sample_outdir.mkdir(parents=True, exist_ok=True)

        # Run per-sample workflow
        sample_command = [
            "cwltool",
            "--cachedir", str(cache_dir),
            "--outdir", str(sample_outdir),
            "--tmpdir-prefix", str(external_tmp) + "/",
            str(per_sample_workflow),
            str(sample_config_file)
        ]

        print(f"   Running fastp + RSEM...")
        sample_result = subprocess.run(sample_command, capture_output=True, text=True)

        if sample_result.returncode != 0:
            print(f"   ❌ Sample {sample_name} failed!")
            if sample_result.stderr:
                # Print last few lines of error
                err_lines = sample_result.stderr.strip().split('\n')
                for line in err_lines[-5:]:
                    print(f"      {line}")
            failed_samples.append(sample_name)
            continue

        # Move outputs to final directories
        print(f"   Moving outputs to final directory...")

        # Find and move files
        for f in sample_outdir.rglob('*'):
            if not f.is_file():
                continue

            fname = f.name

            # Fastp outputs
            if fname.endswith('_R1_clean.fastq.gz') or fname.endswith('_R2_clean.fastq.gz'):
                copy_file_to_dir(f, dir_fastp)
            elif fname.endswith('.html') and 'fastp' in str(f):
                copy_file_to_dir(f, dir_fastp)
            elif fname.endswith('.json') and 'fastp' in str(f):
                copy_file_to_dir(f, dir_fastp)

            # RSEM outputs
            elif fname.endswith('.genes.results'):
                dest = copy_file_to_dir(f, dir_rsem)
                if dest:
                    genes_results_files.append(str(dest))
            elif fname.endswith('.isoforms.results'):
                copy_file_to_dir(f, dir_rsem)

            # BAM files
            elif fname.endswith('.sorted.bam'):
                copy_file_to_dir(f, dir_bams)
            elif fname.endswith('.sorted.bam.bai'):
                copy_file_to_dir(f, dir_bams)

            # Log files
            elif fname.endswith('.rsem.log') or fname.endswith('.rsem.err'):
                copy_file_to_dir(f, dir_logs)
            elif 'Log.' in fname:
                copy_file_to_dir(f, dir_logs)

        # Clean up cache for this sample
        print(f"   Cleaning up cache...")
        deleted = cleanup_sample_cache(cache_dir, sample_name)
        print(f"      Deleted {deleted} cache entries")

        # Clean up sample output directory
        shutil.rmtree(sample_outdir, ignore_errors=True)

        sample_elapsed = time.time() - sample_start
        print(f"   ✅ {sample_name} complete ({sample_elapsed/60:.1f} min)")

    # ========== STEP 3: Run DESeq2 ==========
    if len(genes_results_files) > 0:
        print("\n" + "=" * 70)
        print("Step 3: Running DESeq2 Analysis")
        print("=" * 70)

        # Filter to only successfully processed samples
        successful_samples = [s for s in sample_names if s not in failed_samples]
        successful_groups = [sample_groups[i] for i, s in enumerate(sample_names) if s not in failed_samples]

        deseq_workflow = script_dir / "cwl/tools/deseq-analysis-rsem-standalone.cwl"
        deseq_config = {
            'r_script': {'class': 'File', 'path': str(script_dir / "scripts/deseq_analysis_rsem.R")},
            'gtf_annotation': config_data.get('gtf_annotation'),
            'genes_results_files': [{'class': 'File', 'path': p} for p in genes_results_files],
            'sample_names': successful_samples,
            'sample_groups': successful_groups,
            'output_prefix': config_data.get('run_name', 'deseq2_rsem')
        }

        # Add optional inputs
        if config_data.get('sample_csv_file'):
            deseq_config['sample_csv_file'] = config_data.get('sample_csv_file')
        if config_data.get('geneid2name_file'):
            deseq_config['geneid2name_file'] = config_data.get('geneid2name_file')

        deseq_config_file = external_tmp / "deseq_config.yml"
        with open(deseq_config_file, 'w') as f:
            yaml.dump(deseq_config, f)

        deseq_outdir = external_tmp / "deseq_output"
        deseq_outdir.mkdir(parents=True, exist_ok=True)

        deseq_command = [
            "cwltool",
            "--outdir", str(deseq_outdir),
            "--tmpdir-prefix", str(external_tmp) + "/",
            str(deseq_workflow),
            str(deseq_config_file)
        ]

        print(f"Running DESeq2 on {len(successful_samples)} samples...")
        deseq_result = subprocess.run(deseq_command, capture_output=True, text=True)

        if deseq_result.returncode == 0:
            print("✅ DESeq2 analysis complete!")

            # Move DESeq2 outputs
            print("   Moving DESeq2 outputs...")
            for f in deseq_outdir.rglob('*'):
                if not f.is_file():
                    continue
                fname = f.name
                if fname.endswith('.tsv'):
                    copy_file_to_dir(f, dir_counts, verbose=True)
                elif fname.endswith('.csv') or fname.endswith('.pdf'):
                    copy_file_to_dir(f, dir_analysis, verbose=True)
        else:
            print("❌ DESeq2 analysis failed!")
            if deseq_result.stderr:
                err_lines = deseq_result.stderr.strip().split('\n')
                for line in err_lines[-10:]:
                    print(f"   {line}")

    # ========== Final Summary ==========
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_minutes = elapsed_time / 60
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)

    print("\n" + "=" * 70)
    if len(failed_samples) == 0:
        print("✅ Pipeline completed successfully!")
    else:
        print(f"⚠️  Pipeline completed with {len(failed_samples)} failed samples:")
        for s in failed_samples:
            print(f"   - {s}")

    print(f"\n📁 Results in: {run_directory}/")
    print(f"   ├── {dir_fastp.name}/          (cleaned FASTQ files)")
    print(f"   ├── {dir_rsem.name}/            (gene quantifications)")
    print(f"   ├── {dir_bams.name}/              (sorted BAM files)")
    print(f"   ├── {dir_logs.name}/              (RSEM & STAR log files)")
    print(f"   ├── {dir_analysis.name}/  (DESeq2 results & plots)")
    print(f"   └── {dir_counts.name}/       (raw & normalized counts)")

    print(f"\n⏱️  Total runtime: {elapsed_minutes:.2f} minutes ({hours:02d}:{minutes:02d}:{seconds:02d})")
    print(f"   Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    sys.exit(0 if len(failed_samples) == 0 else 1)


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        prog='long',
        description='RNA-seq Analysis Pipeline - Salmon and RSEM workflows',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                     # Run RSEM pipeline (default)
  %(prog)s rsem                                # Run RSEM pipeline
  %(prog)s salmon                              # Run Salmon pipeline
  %(prog)s rsem --config custom.yml            # Use custom config file
  %(prog)s salmon --path path.yml --config config.yml
        """
    )

    # Add global options that work without subcommand (for default RSEM)
    parser.add_argument(
        '--path',
        type=str,
        help='Path configuration file (default: path_rsem.yml for rsem, path_salmon.yml for salmon)'
    )
    parser.add_argument(
        '--config',
        type=str,
        help='Pipeline configuration file (default: config_rsem.yml for rsem, config_salmon.yml for salmon)'
    )

    subparsers = parser.add_subparsers(dest='command', help='Pipeline to run (default: rsem)')

    # RSEM subcommand (default)
    rsem_parser = subparsers.add_parser(
        'rsem',
        help='Run RSEM + STAR RNA-seq pipeline (default)'
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

    # Salmon subcommand
    salmon_parser = subparsers.add_parser(
        'salmon',
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

    # Parse arguments
    args = parser.parse_args()

    # If no command specified, default to RSEM
    if not args.command:
        run_rsem(args)
    else:
        args.func(args)


if __name__ == '__main__':
    main()
