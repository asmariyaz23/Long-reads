# Long-reads
usage: minion_pipeline.py [-h] --input_path FASTQ_PATH --save_path OUTPUT_DIR
                          [--find_structural_variants] [--resources RESOURCES]
                          [--run_assembly]

Basecalls fast5 files, improves consensus accuracy of an assembly of ONT sequencing reads, aligns to references

required arguments:
  --input_path FASTQ_PATH
                        Fastq files directory paths, comma seperated
  --save_path OUTPUT_DIR
                        Output directory path
  --find_structural_variants
                        Align reads, call structural variants, and annotate
  --resources RESOURCES
                        Resource location
  --run_assembly        Include argument to run assembly program
