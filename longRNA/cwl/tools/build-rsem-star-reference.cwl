#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

doc: |
  Build RSEM + STAR reference with custom genomeSAindexNbases
  Uses a bash script to run both commands in sequence
  
  Save as: cwl/tools/build-rsem-star-reference.cwl

baseCommand: [bash, -c]

inputs:
  genome_fasta:
    type: File
    doc: "Genome FASTA file"
  
  gtf_annotation:
    type: File
    doc: "GTF annotation file"
  
  reference_name:
    type: string
    doc: "Name for the reference"

  sjdboverhang:
    type: int
    doc: "STAR sjdbOverhang (read_length - 1)"

  genomeSAindexNbases:
    type: int
    doc: "STAR genomeSAindexNbases (12 for C. elegans)"
  
  threads:
    type: int?
    doc: "Number of threads"

arguments:
  - valueFrom: |
      ${
        var cmd = "set -e\n";
        cmd += "echo 'Building RSEM transcriptome reference...'\n";
        cmd += "rsem-prepare-reference";
        cmd += " --gtf " + inputs.gtf_annotation.path;
        cmd += " --num-threads " + inputs.threads;
        cmd += " " + inputs.genome_fasta.path;
        cmd += " " + inputs.reference_name;
        cmd += "\n\n";
        cmd += "echo 'Building STAR genome index...'\n";
        cmd += "STAR";
        cmd += " --runMode genomeGenerate";
        cmd += " --runThreadN " + inputs.threads;
        cmd += " --genomeDir .";
        cmd += " --genomeFastaFiles " + inputs.genome_fasta.path;
        cmd += " --sjdbGTFfile " + inputs.gtf_annotation.path;
        cmd += " --sjdbOverhang " + inputs.sjdboverhang;
        cmd += " --genomeSAindexNbases " + inputs.genomeSAindexNbases;
        cmd += "\n\n";
        cmd += "echo 'Reference building complete!'\n";
        return cmd;
      }

outputs:
  reference_dir:
    type: Directory
    outputBinding:
      glob: "."
    doc: "Directory with both RSEM and STAR indices"

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    coresMax: $(inputs.threads)
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.genome_fasta)
      - $(inputs.gtf_annotation)

stdout: combined_reference.log
stderr: combined_reference.err
