#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import argparse
import datetime
import os,sys
import ntpath
import json
import glob
import shutil
import datetime

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "PyLib"]))
from Pipeliner import Pipeliner, Command, run_cmd

import logging
FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
datetime_object = datetime.datetime.now()
datetime_o=str(datetime_object).replace(" ","_")
logging.basicConfig(filename='Minion_pipeline.'+datetime_o+'.log', format=FORMAT, filemode='w', level=logging.DEBUG)
# add a new Handler to print all INFO and above messages to stdout
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
logger.addHandler(ch)

BASEDIR=os.path.dirname(__file__)
UTILDIR=os.path.join(BASEDIR,"util")
GUPPYDIR=os.path.join(UTILDIR,"guppy")
NANOPOLISHDIR=os.path.join(UTILDIR,"nanopolish")
CANU=os.path.join(UTILDIR,"canu")
MINIMAPDIR=os.path.join(UTILDIR,"minimap2")
MERGESCRIPT=os.path.join(UTILDIR,"merge_fastqs.py")
PORECHOP=os.path.join(UTILDIR,"porechop")
HS_BLASTN=os.path.join(UTILDIR,"hs-blastn")
#NANOVAR=os.path.join(UTILDIR,"nanovar")
NGMLR=os.path.join(UTILDIR,"ngmlr")
NPINV=os.path.join(UTILDIR,"npInv1.26.jar")
SNIFFLES=os.path.join(UTILDIR,"sniffles")
ANNOVAR=os.path.join(UTILDIR,"table_annovar.pl")
AnnotSV=os.path.join(UTILDIR,"AnnotSV.tcl")
SMARTDENOVO=os.path.join(UTILDIR,"smartdenovo.pl")
CUTESV=os.path.join(UTILDIR,"cuteSV","cuteSV")
CONVERT_FASTQ_FASTA=os.path.join(UTILDIR,"convert_fastq_fasta.py")
OBINARY_FILE="hg38.fa.counts.obinary"
FASTA_FILE="hg38.fa"
ANNOT_DB="annovar_humandb"
medaka_model="r941_min_fast_g303"
MEDAKA_VARIANT=os.path.join(UTILDIR,"medaka_variant")
canu_spec=os.path.join(UTILDIR,"canu_spec")
minimap2=os.path.join(MINIMAPDIR,"minimap2")
MMI_FILE="hg38.mmi"

class MinIonPipeline:

    def run(self):

        arg_parser = argparse.ArgumentParser(

            description = "Basecalls fast5 files, improves consensus accuracy of an assembly of ONT sequencing reads, aligns to references",
            formatter_class=argparse.RawTextHelpFormatter
            )

        arg_parser._action_groups.pop()
        required = arg_parser.add_argument_group('required arguments')
        required.add_argument("--input_path", dest="fastq_path", type=str, default="", required=True,
                               help="Fastq files directory paths, comma seperated")
        required.add_argument("--save_path", dest="output_dir", type=str, default="", required=True,
                                help="Output directory path")
        '''
        required.add_argument("--guppy_config", dest="guppy_config", type=str, default="",
                                required=True, help=str("Config file to be used for basecalling."+
                                                         " It is located in guppy software directory"+
                                                         " under data/"))
        
        required.add_argument("--run_filt", dest="filt_run", action='store_true',
                              help=str("Include argument to run nanofilt"+
                                       "Skip if reads are already filtered"))
        '''
        required.add_argument("--find_structural_variants", dest="find_structural_variants", 
                               action='store_true',
                              help=str("Align reads, call structural variants, "+
                                       "and annotate"))
        required.add_argument("--resources", dest="resources",
                              help=str("Resource location"))
        required.add_argument("--run_assembly", dest="run_assembly", action='store_true',
                              help=str("Include argument to run assembly "+
                                       "program"))
        args_parsed = arg_parser.parse_args()

        os.system("which python")
        if not os.path.exists(args_parsed.output_dir):
            os.makedirs(args_parsed.output_dir)

        merged_path=os.path.join(args_parsed.output_dir,"merged_reads")
        if not os.path.exists(merged_path):
            os.makedirs(merged_path)
        merged_fastq=os.path.join(args_parsed.output_dir,
                                  "merged_reads",
                                  "merged.fastq")
        compressed_merged_fastq=os.path.join(merged_fastq+".gz")

        filtered_path=os.path.join(args_parsed.output_dir,"filtered_reads")
        if not os.path.exists(filtered_path):
            os.makedirs(filtered_path)
        filtered_fastq=os.path.join(filtered_path,"filtered_reads.fastq")
        filtered_log=os.path.join(filtered_path,"filtered_reads.fastq.log")

        fastq_basecalled_dir=args_parsed.fastq_path

        seq_summary_file=os.path.join(args_parsed.fastq_path,"sequencing_summary.txt")
        stats_path=os.path.join(args_parsed.output_dir,"statistics")
        if not os.path.exists(stats_path):
            os.makedirs(stats_path)
        stats_output=os.path.join(stats_path,"statistics.txt")

        trimmed_path=os.path.join(args_parsed.output_dir,"trimmed_reads")
        if not os.path.exists(trimmed_path):
            os.makedirs(trimmed_path)
        trimmed_output=os.path.join(trimmed_path,"trimmed.fastq")

        if args_parsed.find_structural_variants: 
            nanovar_path=os.path.join(args_parsed.output_dir,"nanovar")
            if not os.path.exists(nanovar_path):
                os.makedirs(nanovar_path)
    
            GENOME_FA=os.path.join(args_parsed.resources,FASTA_FILE)
            GENOME_OBINARY=os.path.join(args_parsed.resources,OBINARY_FILE)
            REF_MMI=os.path.join(args_parsed.resources,MMI_FILE)
    
            nanovar_sv=os.path.join(nanovar_path,"trimmed-hg38.pass.nanovar.vcf")
    
            annot_db=os.path.join(args_parsed.resources,ANNOT_DB)
            annot_path=os.path.join(args_parsed.output_dir,"annotation")
            if not os.path.exists(annot_path):
                os.makedirs(annot_path)
            annotated_vcf=os.path.join(annot_path,"annotation.hg38.vcf")
    
            nanovar_annotsv_path=os.path.join(args_parsed.output_dir,"Annot_SV_NanoVar")
            if not os.path.exists(nanovar_annotsv_path):
                os.makedirs(nanovar_annotsv_path)
            nanovar_annotsv_results=os.path.join(nanovar_annotsv_path,
                                                 "trimmed-hg38.pass.nanovar.annotsv.csv")
    
            align2_path=os.path.join(args_parsed.output_dir,"ngmlr")
            if not os.path.exists(align2_path):
                os.makedirs(align2_path)
            align2_out=os.path.join(align2_path,"trimmed.ngmlr.sam")
    
            sorted_prefix=os.path.join(align2_path,"trimmed.sorted")
    
            sniffles_path=os.path.join(args_parsed.output_dir,"sniffles")
            if not os.path.exists(sniffles_path):
                os.makedirs(sniffles_path)
            sniffles_output=os.path.join(sniffles_path,
                                         "trimmed.sniffles.vcf")
    
            sniffles_annotsv_path=os.path.join(args_parsed.output_dir,"Annot_SV_Sniffles")
            if not os.path.exists(sniffles_annotsv_path):
                os.makedirs(sniffles_annotsv_path)
            sniffles_annotsv_results=os.path.join(sniffles_annotsv_path,
                                                 "trimmed.sniffles.csv")

            cuteSV_path=os.path.join(args_parsed.output_dir,"cuteSV")
            if not os.path.exists(cuteSV_path):
                os.makedirs(cuteSV_path)
            cuteSV_output=os.path.join(cuteSV_path,
                                         "trimmed.cuteSV.vcf")
    
            cuteSV_annotsv_path=os.path.join(args_parsed.output_dir,"Annot_SV_cuteSV")
            if not os.path.exists(cuteSV_annotsv_path):
                os.makedirs(cuteSV_annotsv_path)
            cuteSV_annotsv_results=os.path.join(cuteSV_annotsv_path,
                                                 "trimmed.cuteSV.csv")

        if args_parsed.run_assembly:
            assembly_outdir=os.path.join(args_parsed.output_dir,"canu")
            if not os.path.exists(assembly_outdir):
                os.makedirs(assembly_outdir)

            assembly2_outdir=os.path.join(args_parsed.output_dir,"smartdenovo")
            if not os.path.exists(assembly2_outdir):
                os.makedirs(assembly2_outdir)
            merged_fasta_output=os.path.join(assembly2_outdir,"merged.fasta")
            assembly2_prefix=os.path.join(assembly2_outdir,
                                          "smartdenovo_assembly")

            align2_assembly=os.path.join(args_parsed.output_dir,
                                         "minimap2_aligned_smartdenovo")
            if not os.path.exists(align2_assembly):
                os.makedirs(align2_assembly)
            minimap2_aligned_smartdenovo=os.path.join(align2_assembly,
                                           "minimap2_smartdenovo.sam")
            minimap2_aligned_smartdenovo_bam=os.path.join(align2_assembly,
                                           "minimap2_smartdenovo.bam")

            medaka_variant_smartdenovo=os.path.join(args_parsed.output_dir,
                                         "medaka_variant_smartdenovo")
            if not os.path.exists(medaka_variant_smartdenovo):
                os.makedirs(medaka_variant_smartdenovo)
            medaka_variant_smartdenovo_vcf=os.path.join(medaka_variant_smartdenovo,
                                                 "medaka_variant_smartdenovo.vcf")

            medaka_smartdenovo_annotsv_path=os.path.join(args_parsed.output_dir,
                                                         "Annot_SV_NanoVar")
            if not os.path.exists(medaka_smartdenovo_annotsv_path):
                os.makedirs(medaka_smartdenovo_annotsv_path)
            medaka_smartdenovo_annotsv_output=os.path.join(medaka_smartdenovo_annotsv_path,
                                           "medaka_variant_smartdenovo.vcf.annotated.csv")

        checkpoints_dir = args_parsed.output_dir + "/chckpts_dir"
        checkpoints_dir = os.path.abspath(checkpoints_dir)
        if not os.path.exists(checkpoints_dir):
            os.makedirs(checkpoints_dir)

        ## Construct pipeline
        pipeliner = Pipeliner(checkpoints_dir)

        '''
        ## Basecalling
        cmdstr=str(os.sep.join([GUPPYDIR, "bin","guppy_basecaller"]) +
                      " --input_path " + args_parsed.raw_path +
                      " --save_path " + args_parsed.output_dir +
                      " --config " +  args_parsed.guppy_config +
                      " --recursive " +
                      " --num_callers 1 "+
                      " --cpu_threads_per_caller 1 "+
                      " --qscore_filtering 7")

        pipeliner.add_commands([Command(cmdstr,"basecalling.ok")])
        '''

        ## Merge fastqs
        cmdstr=str("python " + MERGESCRIPT +
                   " --input_dir " + fastq_basecalled_dir +
                   " --merged_output " + merged_fastq)
        pipeliner.add_commands([Command(cmdstr,"merged_fastq.ok")])

        ## Filtering for qscore. Guppy already does this,
        ## however for first few runs, qscore filtering was not
        ## enabled. Hence this needs to be turned on.
        ## Runs only on Python3
        ##if args_parsed.filt_run:
        cmdstr=str("NanoFilt --logfile " + filtered_log+
                   " -q 7 " + merged_fastq +
                   " > " + filtered_fastq)
        pipeliner.add_commands([Command(cmdstr,"filtered.ok")])

        ## Generate statistics. runs only python3
        cmdstr=str("NanoStat --fastq " + compressed_merged_fastq +
                   " --readtype 1D > " + stats_output)
        pipeliner.add_commands([Command(cmdstr,"stats.ok")])

        if args_parsed.find_structural_variants:
            ## Trimming adapters
            cmdstr=str(PORECHOP + " -i " + filtered_fastq +
                       " -o " + trimmed_output)
            pipeliner.add_commands([Command(cmdstr,"trimmed.ok")])
    
            '''
            ## HS_BlastN
            cmdstr=str(HS_BLASTN +" align "+
                       "-reward 2 -penalty -3 "+
                       "-gapopen 0 -gapextend 4 "+
                       "-outfmt 7 " +
                       #"-max_target_seqs 3 "+
                       "-db " + GENOME_FA + " -window_masker_db " +
                       GENOME_OBINARY + " -query " +
                       trimmed_output + " -out " + hs_blastn_fa_out)
            pipeliner.add_commands([Command(cmdstr,"hsblastn.ok")])
            '''
    
    
            ## Alignment2
            cmdstr=str(NGMLR + " -t 4 -r " +
                       GENOME_FA + " -q "+
                       trimmed_output + " -o " + align2_out +
                       " -x ont")
            pipeliner.add_commands([Command(cmdstr,"ngmlr.ok")])
            pipeliner.run()
    
            ## Sort and convert alignment to BAM
            cmdstr=str("samtools sort " + align2_out + " -o "
                        +sorted_prefix+".bam")
            pipeliner.add_commands([Command(cmdstr,"align2.sorted.ok")])
            pipeliner.run()
    
            ## Index BAM file
            cmdstr=str("samtools index " + sorted_prefix+".bam")
            pipeliner.add_commands([Command(cmdstr,"align2.sortedi.indexed.ok")])
            pipeliner.run()
    

            '''
            SV Calling & annotating based pn NGLMR - Sniffles, NVINV, 
            CuteSV
            '''    
            ## SV calling on NGMLR - Sniffles
            cmdstr=str(SNIFFLES + " -m " +
                       sorted_prefix+".bam" + 
                       " --minmapping_qual 5 " +  
                       " -v " + sniffles_output)
            pipeliner.add_commands([Command(cmdstr,"sniffles.ok")])
            pipeliner.run()
    
            ## Annotate Sniffles VCF
            cmdstr=str("$ANNOTSV/bin/AnnotSV/AnnotSV.tcl " + " -SVinputFile " +
                        sniffles_output +
                       " -genomeBuild GRCh38 " +
                       " -outputFile "+ sniffles_annotsv_results)
            pipeliner.add_commands([Command(cmdstr,"sniffles.annotsv.ok")])
            pipeliner.run()

            ## SV calling on NGMLR - CuteSV (meant to be better with short reads)
            cmdstr=str(CUTESV + 
                       " --max_cluster_bias_INS 100 "+
                       " --diff_ratio_merging_INS 0.2 "+
                       " --diff_ratio_filtering_INS 0.6 "+
                       " --diff_ratio_filtering_DEL 0.7 "+
                       " --min_mapq 5 "+
                       sorted_prefix+".bam "+
                       cuteSV_output +" "+ cuteSV_path)
            pipeliner.add_commands([Command(cmdstr,"cuteSV.ok")])
            pipeliner.run()
    
            ## Annotate CuteSV VCF
            cmdstr=str("$ANNOTSV/bin/AnnotSV/AnnotSV.tcl " + " -SVinputFile " +
                       cuteSV_output +
                       " -genomeBuild GRCh38 " +
                       " -outputFile "+ cuteSV_annotsv_results)
            pipeliner.add_commands([Command(cmdstr,"cuteSV.annotsv.ok")])
            pipeliner.run()

            ## Alignment and calling1: Nanovar - Aligns using HS-Blastn and calls SVs
            cmdstr=str("nanovar -t 24 " +
                       trimmed_output + " " + GENOME_FA + " " + nanovar_path)
            pipeliner.add_commands([Command(cmdstr,"nanovar.ok")])
            pipeliner.run()
    
            ## Annotation
            cmdstr=str(ANNOVAR + " " + nanovar_sv + " " +
                       annot_db +
                       " -buildver hg38 "+
                       " -out "+ annotated_vcf + " -remove "+
                       " -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a "+
                       " -operation gx,r,f,f,f -nastring . " +
                       " -vcfinput " +
                       " -polish ")
            pipeliner.add_commands([Command(cmdstr,"annotated.ok")])
            pipeliner.run()
    
            ## Annotate NanoVar VCF
            cmdstr=str("$ANNOTSV/bin/AnnotSV/AnnotSV.tcl " + " -SVinputFile " +
                       nanovar_sv +
                       " -genomeBuild GRCh38 " +
                       " -outputFile "+ nanovar_annotsv_results)
            pipeliner.add_commands([Command(cmdstr,"nanovar.annotsv.ok")])
            pipeliner.run()

        ## Assembly
        if args_parsed.run_assembly:
            '''
            Assemble Using Canu
            '''
            cmdstr=str(CANU +
                       " -p trimmed.canu.assembly " +
                       " -s " + canu_spec +
                       " -d " + assembly_outdir +
                       " genomeSize=3.1g " +
                       " -nanopore-raw " + merged_fastq)
            pipeliner.add_commands([Command(cmdstr,"canu.ok")])
            pipeliner.run()


            ##TODO: work on rest of Canu assembly here



            '''
            Assemble Using Smartdenovo
            '''
            #Convert fastq to fasta
            cmdstr=str("sed -n" + " '1~4s/^@/>/p;2~4p' " +
                       merged_fastq + " > " +
                       merged_fasta_output)
            pipeliner.add_commands([Command(cmdstr,"converted_fastq_fasta.ok")])
            pipeliner.run()

            #Assemble with SMARTDENOVO
            cmdstr=str(SMARTDENOVO + " -p " + assembly2_prefix +
                       " -c 1 -t 4 -k 17 " +
                       merged_fasta_output + " > " + assembly2_prefix +
                       ".mak && "+
                       " make -f " + assembly2_prefix + ".mak")
            pipeliner.add_commands([Command(cmdstr,"smartdenovo.ok")])
            pipeliner.run()

            #Align SMARTDENOVO contigs with Minimap2
            cmdstr=str(minimap2 + " -a " + REF_MMI + " "+ assembly2_prefix+
                       ".fa.gz" +
                       " > " + minimap2_aligned_smartdenovo)
            pipeliner.add_commands([Command(cmdstr,
            "smartdenovo_minimap_aligned.ok")])
            pipeliner.run()

            #Convert Minimap2 sam to bam
            cmdstr=str("samtools sort " + minimap2_aligned_smartdenovo + " -o "
                       + minimap2_aligned_smartdenovo_bam)
            pipeliner.add_commands([Command(cmdstr,
                                    "minimap2_aligned_smartdenovo_bam_converted.ok")])
            pipeliner.run()

            #Index Minimap2 bam
            cmdstr=str("samtools index " + minimap2_aligned_smartdenovo_bam)
            pipeliner.add_commands([Command(cmdstr,
                                    "minimap2_aligned_smartdenovo_bam_converted_indexed.ok")])
            pipeliner.run()

            #Medaka SNP and indel calling on Smartdenovo assembly
            cmdstr=str(MEDAKA_VARIANT +
                       " -m "+ medaka_model +
                       " -f "+ GENOME_FA +
                       " -i "+ minimap2_aligned_smartdenovo_bam +
                       " -o "+ medaka_variant_smartdenovo_vcf)
            pipeliner.add_commands([Command(cmdstr,"medaka_variant_smartdenovo.ok")])
            pipeliner.run()

            #Annotating Medaka variants with AnnotSV
            cmdstr=str("$ANNOTSV/bin/AnnotSV/AnnotSV.tcl " + " -SVinputFile " +
                       medaka_variant_smartdenovo_vcf +
                       " -genomeBuild GRCh38 " +
                       " -outputFile " + medaka_smartdenovo_annotsv_output)
            pipeliner.add_commands([Command(cmdstr,"medaka_smartdenovo_annotsv.ok")])
            pipeliner.run()


if __name__ == "__main__":
    MinIonPipeline().run()

'''    
            npinv_path=os.path.join(args_parsed.output_dir,"npInv")
            if not os.path.exists(npinv_path):
                os.makedirs(npinv_path)
            npinv_output=os.path.join(npinv_path,
                                      "trimmed.npinv.vcf")
    
            npinv_annotsv_path=os.path.join(args_parsed.output_dir,"Annot_SV_npInv")
            if not os.path.exists(npinv_annotsv_path):
                os.makedirs(npinv_annotsv_path)
            npinv_annotsv_results=os.path.join(npinv_annotsv_path,
                                               "trimmed.npinv.csv")
            ## SV calling on NGMLR - npInv
            cmdstr=str("java -jar "+ NPINV +
                       " --output "+ npinv_output +
                       " --input "+ sorted_prefix+".bam")
            pipeliner.add_commands([Command(cmdstr,"npInv.ok")])
            pipeliner.run()
    
            ## Annotate npInv VCF
            cmdstr=str("$ANNOTSV/bin/AnnotSV/AnnotSV.tcl " + " -SVinputFile " +
                       npinv_output +
                       " -genomeBuild GRCh38 " +
                       " -outputFile "+ npinv_annotsv_results)
            pipeliner.add_commands([Command(cmdstr,"npinv.annotsv.ok")])
            pipeliner.run()
'''
