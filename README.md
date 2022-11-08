# Why Unidex?

Unidex is a **uni**versal **de**multiple**x**er written in Python3 that can be applied to a broad range of sequencing data types. Originally built for the considerations of combinatorial indexing, unidex is designed to be flexible in manipulating FASTQ sequence files derived from next-gen sequencing.

A typical workflow generally would start with the sequencing of libraries with an Illumina sequencer, generating the four standard FASTQ files (index1, index2, read1, and read2) from the bcl files using bcl2fastq, and then demultiplexing these reads based on any number of declared index sequence(s) which could be located anywhere within the multiple sequence files.

Unidex is a more flexible tool than other demultiplexers in that it uses a simple language to document and define how to find and return desired sequences from the FASTQ's. The location, amount, and sizes of index sequence used for identifying either specific samples or cells for demultiplexing can vary wildly by technique. Combining these technical indexes and the biological sequence from the multiple read files can also be more complicated than direct concatenations. Users can document their specific 'modes', which capture the specific portions in the FASTQ's comprising index and read sequence, for demultiplexing in a mode file which uses simple language parsable by unidex.


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

Full information about the creation of a mode and mode file is located in a provided example of a mode file in this repository, 'unidex_modes.cfg'. Mode files are tab-separated files with five minimally-required fields. The first is the name of the mode, the second corresponds to the read1 file, the third to the index1 file, the fourth to the index2 file, and the fifth to the read2 file. These fields use the language of 'readN:x', 'indexN:x', and 'null:x' where N is the number designation of the component part and x is the number of bases it covers. The rest of the fields contain the file paths for the declared indexes and non-standard processing flags that can be associated with the mode.

Here is a common example for single-cell combinatorial indexing sequence. The name of this mode is 's3' which corresponds with symmetrical strand sci ('s3') sequencing. The next four fields indicate where the index sequences and read sequence resides within the four files. Here the 10bp index1 is found in the index1 file and corresponds to the i7 index, the 10bp index2 is found in the index2 file and corresponds to the i5 index, and an 8bp index3 is found at the start of read2 which corresponds to the Tn5 sequence. The sequence for read1 comes entirely from the complete read1 file and read2 comes from the remainder of the read2 file skipping the first 8 and 20 bases which correspond to the Tn5 index and adapater sequence respectively.
```
s3	read1:0	index1:10	index2:10	index3:8,null:20,read2:0	index1=/home/groups/oroaklab/src/unidex/indexes/TruSeq_i7_10bp.txt	index2=/home/groups/oroaklab/src/unidex/indexes/Nextera_i5_PCR_10bp.txt	index3=/home/groups/oroaklab/src/unidex/indexes/s3_Tn5_8bp.txt
```

## Running Unidex

```
# Tested example calls.

# Standard calls where the mode, '-M', and run name, '-R', are changing.
# A single holistic mode file, '-m', is being maintained at '/home/groups/oroaklab/demultiplex/unidex_modes.cfg'.
# FASTQ output is standardized in '/home/groups/oroaklab/fastq/' where every run has a directory based on the name of the sequencing run, '-R'.
# Each FASTQ directory has four files named the same: 'Undetermined_S0_I1_001.fastq.gz', 'Undetermined_S0_I2_001.fastq.gz', 'Undetermined_S0_R1_001.fastq.gz', 'Undetermined_S0_R2_001.fastq.gz'.
# Output, '-O', is saved into a directory in '/home/groups/oroaklab/demultiplex/' named after '-R'
python unidex.py -M 220516_iCell8 -R 220520_VH00711_21_AAAMJF5HV \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-r /home/groups/oroaklab/fastq/ \
	-O /home/groups/oroaklab/demultiplex/

python unidex.py -M iCell8_scale36 -R 220603_VH00711_23_AAANL7LHV \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-r /home/groups/oroaklab/fastq/ \
	-O /home/groups/oroaklab/demultiplex/

python unidex.py -M iCell8_scale24 -R 221105_VH00711_42_AACGGTNM5 \
        -m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
        -r /home/groups/oroaklab/fastq/ \
        -O /home/groups/oroaklab/demultiplex/

# example with multiple modes
python unidex.py -M 220916_HTAN_SCR,220916_HTAN_s4 -R 220927_VH00711_37_AACCWJYM5 \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-r /home/groups/oroaklab/fastq/ \
	-O /home/groups/oroaklab/demultiplex/

python unidex.py -M 220916_HTAN_SCR,220930_CEDAR_Human,220930_CEDAR_Mouse -R 220930_VH00711_38_AACCWJFM5 \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-r /home/groups/oroaklab/fastq/ \
	-O /home/groups/oroaklab/demultiplex/

# example with usage of an annotation file
python unidex.py -M sciMET -R 221020_VH00711_40_AACCWJ2M5 \
	-A /home/groups/oroaklab/demultiplex/index_and_annotation_files/quick_annot_make_3_index/annot_files/221020_sciMET.annot \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-r /home/groups/oroaklab/fastq/ \
	-O /home/groups/oroaklab/demultiplex/

# This is an identical call as the one's above but with the FASTQ file names explicitly used instead of guessed.
python unidex.py -M iCell8_scale36 -R 220603_VH00711_23_AAANL7LHV \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg \
	-1 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_R1_001.fastq.gz \
	-2 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_I1_001.fastq.gz \
	-3 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_I2_001.fastq.gz \
	-4 /home/groups/oroaklab/fastq/220603_VH00711_23_AAANL7LHV/Undetermined_S0_R2_001.fastq.gz \
	-O /home/groups/oroaklab/demultiplex/

# A similar call to above but without the '-r' and '-O' flags.
# The full path is in '-R'.
# The prefix name would be taken from that full path as '221105_VH00711_42_AACGGTNM5'.
# The demultiplex output folder would be placed in the same directory as '-R'.
python unidex.py -M iCell8_scale24 -R /home/groups/oroaklab/fastq/221105_VH00711_42_AACGGTNM5 \
	-m /home/groups/oroaklab/demultiplex/unidex_modes.cfg
```


# Contributers and Contact Information

Unidex was originally designed in Perl by Dr. Andrew Adey at Oregon Health and Science University (OHSU). The tool has been translated into Python by James Adler at the Cancer Early Detection Advanced Research Center (CEDAR), a part of the Knight Cancer Research Institue at OHSU. Questions about the tool or help with its usage can be directed to Andrew Nishida, nishidaa@ohsu.edu.

# Known Issues

Relative to the Perl verison, with the use of multiple modes and/or annotation files and using a hamming distance there is a negligible discrepancy of ~1000 reads being differentially assigned in either the Perl version, Python version, or both. In all cases, assigned reads are inflated in the Perl version relative to the Python version and these reads end in the 'Failed' category in the Python version (and not assigned or unassigned due to no annotation). At the least, there is one bug in the original Perl version wherein the total number of output reads surpasses the original number of input read. The following reports demonstrate this bug:
```
/home/groups/oroaklab/demultiplex/220930_VH00711_38_AACCWJFM5/220930_VH00711_38_AACCWJFM5.unidex.report.txt 
/home/groups/oroaklab/demultiplex/220927_VH00711_37_AACCWJYM5/220927_VH00711_37_AACCWJYM5.unidex.report.txt 
```
It is possible there is another bug in the newer Python version but with the Perl version definintely bugged there is no easy way to compare to a known truth. Since this manifests only when there are multiple modes or annotation files and relates directly to hamming distance it is likely the bug in the Perl version relates to the assignment of these alternate barcodes which potentially collide between modes. It is possible the Python version is correct and the Perl version contains the only bugs. Since the output is only minorly changed and these few reads are being put into the 'failed' category and not being falsly assigned downstream to the wrong sample, this discrepency is not being evaluated at the moment.

# Features Coming Soon

- parallelize
- auto-gzip output (but the default python implementation is too slow)
- fold in a 'delayed' or automated mode (maybe move this outside unidex though)
- add additional flags and options for more specific processing -- index-specific hamming, etc.
- add flags and options to allow for different barcode output in the read name
- rewrite read output dictionary to be flexbile in generating any designated number of output read files


# Full Options
```
usage: unidex.py [-h] [-L] [-I REQUEST_MODE_INFO] [-R RUN_FOLDER]
                 [-M MODE_LIST] [-l] [-O OUTPUT_FOLDER]
                 [-d MAX_HAMMING_DISTANCE] [-r FASTQ_FOLDER] -m
                 MODE_CONFIG_FILE [-1 READ1_FILE] [-4 READ2_FILE]
                 [-2 INDEX1_FILE] [-3 INDEX2_FILE] [-A ANNOTATION_FILES]

Demultiplex FASTQs.

optional arguments:
  -h, --help            show this help message and exit

Info options:
  -L, --query_mode_file
                        List modes present in the mode config file. Can
                        specify a different mode file with -m and it will list
                        modes in that file. Can also provide an argument to
                        match to refine list, e.g. 's3'.
  -I REQUEST_MODE_INFO, --request_mode_info REQUEST_MODE_INFO
                        Provide info on one or more comma separated modes as
                        detailed in the specifiedmodes file (-m or default).

Run options:
  -R RUN_FOLDER, --run_folder RUN_FOLDER
                        Run folder where fastq files are present.
  -M MODE_LIST, --mode_list MODE_LIST
                        Mode list - modes must be specified in the modes.cfg
                        file. Modes must be comma separated and will
                        demultiplex in specified order listed.
  -l, --delayed_mode    Delayed mode. Will wait until fastq files are
                        propagated in the specified fastq directory (-r), then
                        will run. Only works when specifying run name, not
                        individual fastq files.

Default options:
  -O OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder (def = Run folder, -R)
  -d MAX_HAMMING_DISTANCE, --max_hamming_distance MAX_HAMMING_DISTANCE
                        Max allowed hamming distance, default 2

Default locations:
  -r FASTQ_FOLDER, --fastq_folder FASTQ_FOLDER
                        Fastq folder full path (def = Run folder, -R.
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
