# Why Unidex?

Unidex is a *uni*versal *de*multiple*x*er written in Python3 that can be applied to a broad range of sequencing data types. Originally built for the considerations of combinatorial indexing, unidex is designed to be flexible in manipulating FASTQ sequence files derived from next-gen sequencing.

A typical workflow generally would start with the sequencing of libraries with an Illumina sequencer, generating the four standard FASTQ files (index1, index2, read1, and read2) from the bcl files using bcl2fastq, and then demultiplexing these reads based on any number of declared index sequence(s) which could be located anywhere within the multiple sequence files.

Unidex is a more flexible tool than other demultiplexers in that it uses a simple language to document and define how to find and return desired sequences from the FASTQ's. The location, amount, and sizes of index sequence used for identifying either specific samples or cells for demultiplexing can vary wildly by technique. Combining these technical indexes and the biological sequence can also be more complicated than direct concatenations. Users can document their specific 'modes' for demultiplexing in a mode file which uses simple language parsable by unidex. Modes capture the specific portions in the FASTQ that comprise index sequence and read sequence.

# Using Unidex

## Installing
Unidex is a stand-alone Python script. It was tested using Python 3.6.8 but should be usable with any version of Python3. The only non-standard Python library that needs to be installed is 'exrex', `pip install exrex`.

If needed, here is a simple way to setup an environment for running unidex.py:
```
# make environment
python3 -m venv unidex_env

# activate environment
source /home/groups/oroaklab/nishida/env/unidex_env/bin/activate
pip install --upgrade pip
pip install --upgrade setuptools
pip install exrex
```

## Setting up a Mode File

Full information about the creation of a mode and mode file is located in a provided example of a mode file in this repository, 'unidex_modes.cfg'. Mode files are tab-separated files with five minimally-required fields. The first is the name of the mode, the second corresponds to the read1 file, the third to the index1 file, the fourth to the index2 file, and the fifth to the read2 file. These fields use the language of 'readN:x', 'indexN:x', and 'null:x' where N is the number designation of the component part and x is the number of bases it covers. The rest of the fields contain the file paths for the declared indexes.

Here is a common example for single-cell combinatorial indexing sequence. The name of this mode is 's3' which corresponds with symmetrical strand sci ('s3') sequencing. The next four fields indicate where the index sequences and read sequence resides within the four files. Here the 10bp index1 is found in the index1 file and corresponds to the i7 index, the 10bp index2 is found in the index2 file and corresponds to the i5 index, and an 8bp index3 is found at the start of read2 which corresponds to the Tn5 sequence. The sequence for read1 comes entirely from the complete read1 file and read2 comes from the remainder of the read2 file skipping the first 8 and 20 bases which correspond to the Tn5 index and adapater sequence respectively.
```
s3	read1:0	index1:10	index2:10	index3:8,null:20,read2:0	index1=/home/groups/oroaklab/src/unidex/indexes/TruSeq_i7_10bp.txt	index2=/home/groups/oroaklab/src/unidex/indexes/Nextera_i5_PCR_10bp.txt	index3=/home/groups/oroaklab/src/unidex/indexes/s3_Tn5_8bp.txt
```

## Running Unidex

```
# placeholder

# example run that reproduces exact counts to perl version
python unidex.py -M 220516_iCell8 -m unidex_modes.cfg -R 220520_VH00711_21_AAAMJF5HV -r /home/groups/oroaklab/fastq/220520_VH00711_21_AAAMJF5HV/ \
	-1 /home/groups/oroaklab/fastq/220520_VH00711_21_AAAMJF5HV/Undetermined_S0_R1_001.fastq.gz \
	-2 /home/groups/oroaklab/fastq/220520_VH00711_21_AAAMJF5HV/Undetermined_S0_I1_001.fastq.gz \
	-3 /home/groups/oroaklab/fastq/220520_VH00711_21_AAAMJF5HV/Undetermined_S0_I2_001.fastq.gz \
	-4 /home/groups/oroaklab/fastq/220520_VH00711_21_AAAMJF5HV/Undetermined_S0_R2_001.fastq.gz \
	-O test_py
	
# example run that reproduces exact counts to perl version
python unidex.py -M iCell8_scale36 -m unidex_modes.cfg -R 220603_VH00711_23_AAANL7LHV -r /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/ \
	-1 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_R1_001.fastq.gz \
	-2 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_I1_001.fastq.gz \
	-3 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_I2_001.fastq.gz \
	-4 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_R2_001.fastq.gz \
	-O test_py
```

# Contributers and Contact Information

Unidex was originally designed in Perl by Dr. Andrew Adey at Oregon Health and Science University (OHSU). The tool has been translated into Python by James Adler at the Cancer Early Detection Advanced Research Center (CEDAR), a part of the Knight Cancer Research Institue at OHSU. Questions about the tool or help with its usage can be directed to Andrew Nishida, nishidaa@ohsu.edu.

# Features Coming Soon
- parallelization
- delayed flags
- xoxo

# Full Options
```
usage: unidex.py [-h] [-L] [-I REQUEST_MODE_INFO] [-R RUN_FOLDER]
                 [-M MODE_LIST] [-l] [-O OUTPUT_FOLDER]
                 [-d MAX_HAMMING_DISTANCE] [-r FASTQ_FOLDER] -m
                 MODE_CONFIG_FILE [-1 READ1_FILE] [-4 READ2_FILE]
                 [-2 INDEX1_FILE] [-3 INDEX2_FILE] [-A ANNOTATION_FILES]

Demultiplex fastq

optional arguments:
  -h, --help            show this help message and exit

Info options:
  -L, --query_mode_file
                        List modes present in the mode config file.Can specify
                        a different mode file with -m and it will list modes
                        in that file.Can also provide an argument to match to
                        refine list, e.g. 's3'
  -I REQUEST_MODE_INFO, --request_mode_info REQUEST_MODE_INFO
                        Provide info on one or more comma separated modes as
                        detailed in the specifiedmodes file (-m or default).

Run options:
  -R RUN_FOLDER, --run_folder RUN_FOLDER
                        Run Folder (where fastq files are present)
  -M MODE_LIST, --mode_list MODE_LIST
                        Mode list - modes must be specified in the modes.cfg
                        file.Modes must be comma separated and will
                        demultiplex in specified order listed.
  -l, --delayed_mode    Delayed mode. Will wait until fastq files are
                        propagatedin the specified fastq directory (-r), then
                        will run.Only works when specifying run name, not
                        individual fastq files.

Default options:
  -O OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder (def = run name, -R)
  -d MAX_HAMMING_DISTANCE, --max_hamming_distance MAX_HAMMING_DISTANCE
                        Max allowed hamming distance

Default locations:
  -r FASTQ_FOLDER, --fastq_folder FASTQ_FOLDER
                        Fastq folder full path.
  -m MODE_CONFIG_FILE, --mode_config_file MODE_CONFIG_FILE
                        Mode config file

Fastq input:
  -1 READ1_FILE, --read1_file READ1_FILE
                        Read 1 fastq
  -4 READ2_FILE, --read2_file READ2_FILE
                        Read 2 fastq
  -2 INDEX1_FILE, --index1_file INDEX1_FILE
                        Index 1 fastq
  -3 INDEX2_FILE, --index2_file INDEX2_FILE
                        Index 2 fastq

Other options:
  -A ANNOTATION_FILES, --annotation_files ANNOTATION_FILES
                        Annotation file(s), comma separated with mode
                        specifiedIf only one mode is specified, then it will
                        default to that mode.[mode1]=[annot_file1],[mode2]=[an
                        not_file2],etc... ORFirst column of annot file
                        designated mode for that annot
```