#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import argparse
import subprocess


def convert_fq_fa(fastq_file,fastq_output):
    cmd=("sed -n" + " '1~4s/^@/>/p;2~4p' " + 
         fastq_file + " > " + fasta_file)
    subprocess.check_call(cmd, shell=True) 

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--fastq_file", help="Fastq file to convert")
    parser.add_argument("--fasta_output", help="Output fasta")
    args = parser.parse_args()



