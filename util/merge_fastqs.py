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
import gzip

def merge_files(f_list,concatenated_file_path):
    with open(concatenated_file_path,'wb') as wfd:
        for fastq in f_list:
            with open(fastq,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

def compress(merged_output):
    compressed_file=merged_output+".gz"
    with open(merged_output, 'rb') as f_in:
        with gzip.open(compressed_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--input_dir", help="Directory path where basecalled fastqs reside")
    parser.add_argument("--merged_output", help="Final concatenated fastq file")
    args = parser.parse_args()

    input_list=(args.input_dir).split(",")
    all_fastq_files=[]
    for il in input_list:
        all_files=os.listdir(il)
        if "pass" in all_files:
            pass_dir=os.path.join(il,"pass")
            passed_f=os.listdir(pass_dir)
            f_list=[os.path.join(pass_dir,f) for f in passed_f if f.endswith(".fastq")]
        else:
            f_list=[os.path.join(il,f) for f in all_files if f.endswith(".fastq")]
        all_fastq_files.extend(f_list)
    merge_files(all_fastq_files,args.merged_output)
    compress(args.merged_output)
