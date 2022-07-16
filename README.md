# Kmer assembler
## Table of Contents

   * [Introduction](#introduction)

   * [Usage](#usage)
      * [Simple Usage](#Usage)
      * [Inputs](#inputs)
      * [Outputs](#outputs)
   * [Citations](#citations)
   * Known issues


### Introduction

Ran_Kmer_assembler is a command line interface for extend short sequences by kmer overlaps. Although it is designed to work at a small scale, but it should work with sequences of any length (Not tested). Ran_Kmer_assembler works as a wrapper on top of the core assembler (kmer_asm2.pl).

One common scenario for this assmbler is after you harvesting a certain sets of kmers (e.g. most abundant kmers in a genome can be used to construct the repeat unit), you want to extend the short kmers into a much longer sequences that could be used for other analysis.

## Why another assembler?
I need a quick tool to assemble short sequences (length<30) into repeat unit. And it is a good practice to understand how the basic assembly is done.

### Usage

`perl Ran_Kmer_assembler.pl SEQ_FILE.fa OVERLAP_SIZE`
Only two parameters are needed. `SEQ_FILE.fa` is a fasta format sequence file, and `OVERLAP_SIZE` is the length to allow a merge. And the overlap_size should be no more than the length of input sequences.

### Output
Each round will be in its own folder. The assembler will stop when there is no sequence to merge. So the results from the last round (n) should be the same as the n-1 round.


### Installation
All you need is Perl installed in your OS

### Known issues and missing parts
There is no extension check outside the kmer overlap, which means this assembler does not support calling alleles and error could occur due to unexpected overlaps.
There is no validation or systematic check on the error rate.
Since the inputs are hashed inorderly, the results would be different each time. I might update the statistic calculation upon further request.
