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

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "PyLib"]))
from Pipeliner import Pipeliner, Command, run_cmd

import logging
FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(filename='Minion_pipeline.log', format=FORMAT, filemode='w', level=logging.DEBUG)
# add a new Handler to print all INFO and above messages to stdout
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
logger.addHandler(ch)

BASEDIR=os.path.dirname(__file__)
UTILDIR=os.path.join(BASEDIR,"util")
GUPPYDIR=os.path.join(UTILDIR,"guppy")
NANOPOLISHDIR=os.path.join(UTILDIR,"nanopolish")
CANUDIR=os.path.join(UTILDIR,"canu")
MINIMAPDIR=os.path.join(UTILDIR,"minimap2")

class MinIonPipeline:

    def run(self):
        
        arg_parser = argparse.ArgumentParser(
            
            description = "Basecalls fast5 files, improves consensus accuracy of an assembly of ONT sequencing reads, aligns to references",
            formatter_class=argparse.RawTextHelpFormatter
            )
        
        arg_parser._action_groups.pop()
        required = arg_parser.add_argument_group('required arguments')
        required.add_argument("--input_path", dest="raw_path", type=str, default="", required=True, 
                               help="Raw fast5 files directory path")
        required.add_argument("--save_path", dest="output_dir", type=str, default="", required=True,
                                help="Output directory path")
        required.add_argument("--guppy_config", dest="guppy_config", type=str, default="", 
                                required=True, help=str("Config file to be used for basecalling."+ 
                                                         " It is located in guppy software directory"+
                                                         " under data/"))
        args_parsed = arg_parser.parse_args()
        
        if not os.path.exists(args_parsed.output_dir):
            os.makedirs(args_parsed.output_dir)

        checkpoints_dir = args_parsed.output_dir + "/chckpts_dir"
        checkpoints_dir = os.path.abspath(checkpoints_dir)
        if not os.path.exists(checkpoints_dir):
            os.makedirs(checkpoints_dir)

        ## Construct pipeline
        pipeliner = Pipeliner(checkpoints_dir)
   
        ## Basecalling
        cmdstr=str(os.sep.join([GUPPYDIR, "bin","guppy_basecaller"]) +
                      " --input_path " + args_parsed.raw_path +
                      " --save_path " + args_parsed.output_dir +
                      " --config " +  args_parsed.guppy_config)
                                   
        pipeliner.add_commands([Command(cmdstr,"basecalling.ok")])
        pipeliner.run()
                
        ## Merge fastqs
        basecall_op=args_parsed.output_dir
	all_files=os.listdir(basecall_op)
        fastqs=[f for f in all_files if f.endswith(".fastq")]

        merged_fastq=os.path.join(args_parsed.output_dir,'guppy.basecalled.fastq')
        if not os.path.exists(merged_fastq):
	    with open(merged_fastq,'wb') as wfd:
		for fastq in fastqs:
                    full_path=os.path.join(args_parsed.output_dir,fastq)
		    with open(full_path,'rb') as fd:
			shutil.copyfileobj(fd, wfd)

        ## Indexing reads 
        if os.sep.join([basecall_op,"sequencing_summary.txt"]):
            if fastqs:
                cmdstr = str( os.sep.join([NANOPOLISHDIR,"nanopolish"]) +
                              " index -d " + args_parsed.raw_path +
                              " " +merged_fastq)

        pipeliner.add_commands([Command(cmdstr,"nanopolish_index.ok")]) 

        ## Compute draft genome
        
        cmdstr = str(os.sep.join([CANUDIR,"canu"]) + 
                     " -p human -d "+ args_parsed.output_dir + 
                     " genomeSize=3.1gb " + #"stopOnLowCoverage=1" +
                     " -nanopore-raw " + merged_fastq)
         
        pipeliner.add_commands([Command(cmdstr,"draft_genome.ok")]) 

        ## Improve the assembly

        assembly=os.path.join(args_parsed.output_dir,
                              "human.contigs.fasta")
        sorted_bam=os.path.join(args_parsed.output_dir,
                                "reads.sorted.bam")
        temp_reads=os.path.join(args_parsed.output_dir,
                                "reads.tmp")
        cmdstr = str(os.sep.join([MINIMAPDIR,"minimap2"])+
                     " -ax map-ont -t 8 "+
                     assembly + " " +
                     merged_fastq +
                     " | samtools sort -o "+ 
                     sorted_bam +
                     " -T "+ temp_reads)
        pipeliner.add_commands([Command(cmdstr,"improve_assem.ok")])

        
        ## Index bam file
        cmdstr = str("samtools index " + sorted_bam)
        pipeliner.add_commands([Command(cmdstr,"bam_index.ok")])

        ## Large files chuncks processed into variants finding algo
        nano_result=os.path.join(args_parsed.output_dir,
                                 "nanopolish.results")
        nano_vcf= os.path.join(args_parsed.output_dir,
                               "polished.vcf")
        cmdstr = str(os.sep.join([NANOPOLISHDIR,"scripts","nanopolish_makerange.py "])
                                 + assembly +" | "+ "parallel --results " 
                                 + nano_result + " -P 8 " 
                                 + os.path.join(NANOPOLISHDIR,"nanopolish")
                                 + " variants --consensus "
                                 + "-o "+nano_vcf+ " -w {1}" 
                                 + " -r " +merged_fastq 
                                 + " -b " +sorted_bam 
                                 + " -g "+assembly+ " -t 4 " 
                                 + "--min-candidate-frequency 0.1")

        pipeliner.add_commands([Command(cmdstr,"variant_calling.ok")])
        
        pipeliner.run()

if __name__ == "__main__":
    MinIonPipeline().run()
