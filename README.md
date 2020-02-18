# ONT Long-reads (SV calling and assembly)
```
usage: minion_pipeline.py [-h] --input_path FASTQ_PATH --save_path OUTPUT_DIR
                          [--find_structural_variants] [--resources RESOURCES]
                          [--run_assembly]

Improves consensus accuracy of an assembly of ONT sequencing reads, aligns to references and calls structural variants

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
 ```
 
This pipeline uses the following 3rd party tools for aligning and SV calling:

1) NanoFilt 
2) NanoStat
3) Porechop
4) NGMLR
5) Sniffles
6) CuteSV
7) NanoVar
8) AnnotSV

This pipeline uses the following 3rd party tools for assembling and aligning:

1) Canu
2) Smartdenovo
3) minimap2
4) Medaka Variant
5) AnnotSV
