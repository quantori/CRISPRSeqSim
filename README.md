### CRISPRSeqSim

CrispTestDataGenerator is a tool that generates variations of the initial and desired sequence
(simulation of CRISPR/Cas editing) and provides a report about the mutations and their positions.
The data produced by CrispTestDataGenerator can be used for testing purposes.

How to use:
1. In the Config file adjust the settings to generate the sequence of desired length with corresponding mutations;

{
      "initLength": 15, # the length of the desired sequence; must be divisible by 3, if "isCoding" parameter is true
      "isCoding": true, # a bool value; if "false", the "initLength" can be of any size
      "mutationSite": [15], # the sites where mutations in the "mutated sequence" (a sequence of interest after CRISPR/Cas editing) are to be made
      "mutationLength": [1], # the type of the mutation in the site. Insertion - positive integer; deletion - negative integer, substitution - 0;
      "direction": [true, true, true], # false if the mutation is directed to the left (in case with deletions)
      "quality": 50.0 # the percentage of "corrupted" mutated sequences in the reads, that do not correspond to the mutated sequence
    }
2. In the "sequences_generator.py" file (method "corrupt_generated") adjust the count. Count is the overall number of the
sequences in the file with reads. You can also adjust the number of the sequences with no mutations at all (initial sequence) in the
c_int variable:

    [def corrupt_generated(self, orig_sequence, sequence):
        count = 100 # can be adjusted
        c_int = random.randint(0, c_count // 2)] # can be adjusted

How it works:
Output

The file with sequences: output.fasta
1. The initial (reference sequence) of the adjusted length is generated randomly;
2. The mutated sequence (with mutations on mutation sites from the Config file) is generated based on the initial sequence;


The file with CRISPR/Cas reads: output.fastq
The file contains the set number of counts. Based on the quality score the certain percentage of the sequences is corrupted
(has other mutations, than the mutated sequence). Also the certain amount of the sequences are unchanged at all (initial or
reference sequence). The number of such sequence can be found in the report file.

sizeScale - the field that specifies the range within which the length of the sequences in the FASTQ file should be modified;
100% - non modified;

scaleDirection: -1 - the length modification is performed to the left side;
0 - from the both sides;
1- to the right side;
