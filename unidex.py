#!/usr/bin/env python

# unidex initial python build

# import modules
import argparse
import sys
import os
import logging
import json
import re
import gzip
import exrex
import itertools


def parse_args():
    parser = argparse.ArgumentParser(description = "Demultiplex fastq")

    # Info options
    parser.add_argument("-L", "--query_mode_file", action = 'store_true', help = "List modes present in the mode config file."
                            "Can specify a different mode file with -m and it will list modes in that file."
                            "Can also provide an argument to match to refine list, e.g. 's3'")
    parser.add_argument("-I", "--request_mode_info", help = "Provide info on one or more comma separated modes as detailed in the specified"
                            "modes file (-m or default).")

    # Run options
    parser.add_argument("-R", "--run_folder", help = "Run Folder (where fastq files are present)")
    parser.add_argument("-M", "--mode_list", type = str, help = "Mode list - modes must be specified in the modes.cfg file."
                            "Modes must be comma separated and will demultiplex in specified order listed.")
    parser.add_argument("-l", "--delayed_mode", action = 'store_true', help = "Delayed mode. Will wait until fastq files are propagated"
                            "in the specified fastq directory (-r), then will run."
                            "Only works when specifying run name, not individual fastq files.")

    # Default options
    parser.add_argument("-O", "--output_folder", help = "Output folder (def = run name, -R)")
    parser.add_argument("-d", "--max_hamming_distance", help = "Max allowed hamming distance", type = int, default = 2)

    # Default locations
    parser.add_argument("-m", "--mode_config_file", help = "Mode config file", required = True)

    # Fastq input (default = auto detect):
    parser.add_argument("-1", "--read1_file", help = "Read 1 fastq")
    parser.add_argument("-4", "--read2_file", help = "Read 2 fastq")
    parser.add_argument("-2", "--index1_file", help = "Index 1 fastq")
    parser.add_argument("-3", "--index2_file", help = "Index 2 fastq")

    # Other options
    parser.add_argument("-A", "--annotation_files", type = str, help = "Annotation file(s), comma separated with mode specified"
                            "If only one mode is specified, then it will default to that mode."
                            "[mode1]=[annot_file1],[mode2]=[annot_file2],etc... OR"
                            "First column of annot file designated mode for that annot")
    
    return parser.parse_args()


# TODO: can add additional funcitonality here to look for valid mode file
def validate_mode_file_exists(mode_file:str) -> bool:
    """
    Exits if mode file doesn't exist. Returns True if path to mode file is valid.

    Parameters:
    -----------
    mode_file : str
        Path to mode config file
    """

    if not os.path.exists(mode_file):
        sys.exit("Please define valid mode file using -m flag.")
    return True

        
def print_available_modes(mode_file:str) -> None:
    """
    Prints available modes to console and exits program.

    Parameters:
    -----------
    mode_file : file
        Path to mode config file
    """
    validate_mode_file_exists(mode_file)
    print("Mode file specified:\n\t{}".format(mode_file))
    print("Modes available:")
    with open(mode_file) as f:
        while True:
            line = f.readline()
            if not line: # break at end of file
                break
            if line.startswith("#"): # skip headers
                continue
            print("\t{}\n".format(line.split('\t')[0].strip()), sep = "", end = "") # print mode
    sys.exit(0) # exit without error
    return None


def print_mode_details(mode_file:str, mode:str) -> None:
    """
    Prints details of user-specified mode and exits program.

    Parameters:
    -----------
    mode_file : str
        Path to mode config file.
    mode : str
        User-specified mode of interest.
    """
    validate_mode_file_exists(mode_file)
    print("Mode file specified:\n\t{}".format(mode_file))
    print("User-specified mode of interest:\n\t{}".format(mode))
    print("Mode details:")
    
    # find specified mode in mode file and print details
    with open(mode_file, "r") as f:
        while True:
            line = f.readline()
            if not line: # break at end of file
                break
            if line.startswith(mode):
                # TODO: should be able to add more flexibility here
                mode_details = line.split("\t")
                print("\tRead1 = {}".format(mode_details[1]))
                print("\tIndex1 = {}".format(mode_details[2]))
                print("\tIndex2 = {}".format(mode_details[3]))
                print("\tRead2 = {}".format(mode_details[4]))
                print("\tIndex Files:")
                for index_file in mode_details[5:]:
                    # TODO: add validation of file existence
                    print("\t\t{}".format(index_file.strip()))
                sys.exit(0)
    # if the function makes it here the mode doesn't exist
    sys.exit("Mode {} does not exist in mode config file {}".format(mode, mode_file))
    return None


def validate_run_options(args) -> None:
    """
    Validates combination fo user-specified run options. Exits program if invalid.

    Parameters:
    -----------
    args : type
        Arguments specified by user
    """
    if args.run_folder is None:
        if args.read1_file is not None or args.read2_file is not None or args.index1_file is not None or args.index2_file is not None:
            if args.output_folder is not None:
                sys.exit("ERROR: When not supplying -R and instead specifying individual fastq files, an output name must be provided (-O).")
        else:
            sys.exit("ERROR: If -R is not specified, each input fastq must be specified (min of 1)!")
    if args.mode_list is None:
        sys.exit("ERROR: Modes list must be specified with one or more modes!")
    return None


def parse_comma_separated_inputs(comma_separated_input_string:str) -> list:
    """
    Parses the comma separated string of mode or annotation files passes as input

    Parameters:
    -----------
    comma_separated_input_string : str
        A string of inputs separated by commas - modes or annot files

    Returns:
    --------
    list
        List containing each specified mode
    """
    # split mode string on commas
    return comma_separated_input_string.split(',')


def validate_mode_count(mode_list:list, mode_dict:dict, mode_config_file:str):
    """
    Validates that all modes were added to the dictionary.
    """
    if not all(mode in mode_dict for mode in mode_list):
        # TODO:
        sys.exit("Dictionary does not contain all specified modes:\n"
                    "mode_list: {}\n"
                    "mode_dict: {}\n"
                    "ERROR: Please validate modes are specified in config file: {}".format(mode_list, mode_dict.keys(), os.path.realpath(mode_config_file)))


# TODO: define the aspects of the dictionary in the docstrings
def generate_mode_dict(mode_list:list, mode_config_file:str) -> dict:
    """
    Processes specified modes for current run and stores mode characteristics in dictionary.

    Parameters:
    -----------
    mode_list : list
        List of modes specified by user
    mode_config_file : str
        Path to mode configuration file specified by user

    Returns:
    --------
    mode_dict : dict
        Dictionary containing the characterisits of each specified mode.
    """

    logging.info("Generating mode dictionary")

    mode_dict:dict = {} # instantiate mode dictionary
    modes_added:int = 0 # instantiate mode tracker

    # read through config file and store entries specified in modes list
    open_mode_config_file = open(mode_config_file, 'r')

    while modes_added < len(mode_list):
        line = open_mode_config_file.readline()
        if not line: # break at end of file
            break
        if line.startswith("#"): # skip header / description lines
            continue
        mode_components = line.strip().split() # parse mode components / characteristics
        
        # extract mode components if in specified mode list
        if mode_components[0] in mode_list:
            mode_dict = add_mode_to_dict(mode_components, mode_dict)
            modes_added += 1

    open_mode_config_file.close() # close mode config

    validate_mode_count(mode_list, mode_dict, mode_config_file)

    logging.info("Mode dicitionary successfully created") # log dictionary
    for mode in mode_dict:
        logging.info("\t{}:".format(mode))
        for component in mode_dict[mode]:
            logging.info("\t\t{}: {}".format(component, mode_dict[mode][component]))

    return mode_dict


# TODO: may need to come back and address additional mode variations - this works for the s3 mode but may fail for the 10x mode
# TODO: can make much more flexible
def add_mode_to_dict(mode_components:list, mode_dict:dict) -> dict:
    """
    Adds new mode to mode dict

    Parameters:
    -----------
    mode_components : list
        Characteristics of the mode to be added to the mode dictionary

    mode_dict : dict
        Dictionary containing all user specified modes for current run

    Returns:
    --------
    mode_dict : dict
        Returns the mode dictionary with the new mode addition
    """
  
    # extract components of mode
    mode, read1, index1, index2, read2, index_files = mode_components[0], mode_components[1], mode_components[2], mode_components[3], mode_components[4], mode_components[5:]

    # instantiate mode and compoenents within mode dictionary
    mode_dict[mode] = {}
    mode_dict[mode]['read1'] = {part.split(":")[0]: int(part.split(":")[1]) for part in read1.split(",")}
    mode_dict[mode]['index1'] = {part.split(":")[0]: int(part.split(":")[1]) for part in index1.split(",")}
    mode_dict[mode]['index2'] = {part.split(":")[0]: int(part.split(":")[1]) for part in index2.split(",")}
    mode_dict[mode]['read2'] = {part.split(":")[0]: int(part.split(":")[1]) for part in read2.split(",")}

    # add index file paths to dictionary - ignoring special flags
    mode_dict[mode]['index_file_paths'] = {part.split("=")[0]: part.split("=")[1] for part in index_files if not re.search('index[0-9]_', part.split("=")[0])}

    return mode_dict


def generate_annotation_dict(annotation_files_list:list, mode_list:list=None) -> dict:
    """
    Process annotation file(s).

    Parameters:
    -----------
    annotation_files_list : list
        List of all annotation files with specified mode. Format is [mode1]=[annot_file1]
    mode_list : list, default = None
        List of corresponding mode designation for each annotation file. Default value is None, as mode could be specified in annotation file.
    
    Returns:
    --------
    annotation_dict : dict
    """
    # instantiate empty annotation dictionary
    annotation_dict:dict = {}

    # process each annotation file
    for i, annotation_file in enumerate(annotation_files_list):
        logging.info("Processing annotation file: {}".format(annotation_file))
        
        mode_in_annot_file:bool = False # default is for mode to specified at command line
        lines_processed:int = 0 # instantiate line tracker

        # mode specified in command line
        if "=" in annotation_file:
            mode, annotation_file = annotation_file.split("=")
        elif mode_list is not None:
            mode = mode_list[i]

        open_annotation_file = open(annotation_file, 'r')
        # process first line to assess if mode is specified in annot file
        line = open_annotation_file.readline().strip().split()
        if len(line) == 3:
            logging.info("Modes specified in annot file: {}".format(annotation_file))
            mode_in_annot_file = True
            mode, cellID, annot = line
        else:
            cellID, annot = line

        # add mode to annotation dict
        if mode not in annotation_dict:
            annotation_dict[mode] = {}
        annotation_dict[mode][cellID] = annot

        # process entire file
        while True:
            line = open_annotation_file.readline().strip().split()            
            if not line: # break if end of file
                break
            if mode_in_annot_file: # may be multiple modes within annotation file
                mode, cellID, annot = line # extract all three columns
                if mode not in annotation_dict:
                    annotation_dict[mode] = {}
            else:
                cellID, annot = line # extract only two columns if mode specified in command line        
            annotation_dict[mode][cellID] = annot # add instance to dictionary
        
        logging.info("Total lines processed for annotation file '{}': {}".format(annotation_file, len(annotation_dict[mode])))
        open_annotation_file.close()

    return annotation_dict


# TODO: build this out
def execute_delayed_mode():
    sys.exit("Delayed mode not currently implemented.")


def define_input_files(args) -> tuple:
    """
    Defines read and index files. 

    Parameters:
    -----------
    args : type
        All user-defined arguments, including optionally specified read and index file paths.

    Returns:
    --------
    read1_file : str
        R1 fastq file
    read2_file : str
        R2 fastq file
    index1_file : str
        I1 fastq file
    index2_file : str
        I2 fastq file
    """
    # instantiate standard files
    read1_file:str = ""
    read2_file:str = ""
    idnex1_file:str = ""
    index2_file:str = ""
    
    # define read 1
    if args.read1_file is not None:
        read1_file = os.path.abspath(args.read1_file)
    elif args.run_folder is not None:
        read1_file = os.path.abspath(os.path.join(args.run_folder, "Undetermined_S0_L001_R1_001.fastq.gz"))
        if not os.path.exists(read1_file):
            sys.exit("ERROR: Read1 file does not exist {}".format(read1_file))
    else:
        sys.exit("ERROR: Either read1 file (-1 flag) OR run folder (-R flag) must be defined.")

    # define read 2
    if args.read2_file is not None:
        read2_file = os.path.abspath(args.read2_file)
    elif args.run_folder is not None:
        read2_file = os.path.abspath(os.path.join(args.run_folder, "Undetermined_S0_L001_R2_001.fastq.gz"))
        if not os.path.exists(read2_file):
            logging.info("No read 2 file detected at path {}.\nMoving forward with single end read.".format(read2_file))
            read2_file = ""

    # define index 1
    if args.index1_file is not None:
        index1_file = os.path.abspath(args.index1_file)
    elif args.run_folder is not None:
        index1_file = os.path.abspath(os.path.join(args.run_folder, "Undetermined_S0_L001_I1_001.fastq.gz"))
        if not os.path.exists(index1_file):
            logging.info("No index 1 file detected at path {}".format(index1_file))
            index1_file = ""

    # define index 2    
    if args.index2_file is not None:
        index2_file = os.path.abspath(args.index2_file)
    elif args.run_folder is not None:
        index2_file = os.path.abspath(os.path.join(args.run_folder, "Undetermined_S0_L001_I2_001.fastq.gz"))
        if not os.path.exists(index2_file):
            logging.info("No index 2 file detected at path {}".format(index2_file))
            index2_file = ""

    # print existing files to log
    logging.info("Read 1 file: {}".format(read1_file if read1_file else "Not defined or doesn't exist"))
    logging.info("Read 2 file: {}".format(read2_file if read2_file else "Not defined or doesn't exist"))
    logging.info("Index 1 file: {}".format(index1_file if index1_file else "Not defined or doesn't exist"))
    logging.info("Index 2 file: {}".format(index2_file if index2_file else "Not defined or doesn't exist"))

    return read1_file, read2_file, index1_file, index2_file


# TODO: define dictionary components within docstrings
def generate_expected_index_dict(mode_dict:dict, hamming_distance:int) -> dict:
    """
    Generates expected index dictionary using expected index files for each mode in config file.

    Parameters:
    -----------
    mode_dict : dict
        Mode dictionary containing modes as keys and mode characateristics as values.
    hamming_dist : 2
        Maximum allowable hamming distance for each index.

    Returns:
    --------
    expected_index_dict : dict
        Dictionary containing modes as keys with values as subdictionaries containing index name
        as key and a list of expected indexes as values.
    """

    # instantiate dict
    expected_index_dict:dict = {}

    # extract modes and indexes from mode_dict and defined index files
    for mode in mode_dict:
        expected_index_dict[mode] = {}
        for index in mode_dict[mode]['index_file_paths']:
            with open(mode_dict[mode]['index_file_paths'][index], "r") as f:
                expected_index_dict[mode][index] = {line.split('\t')[-1]:line.split('\t')[-1] for line in f.read().split('\n') if line}
    
    # generate hamming distance possibilities
    for index in expected_index_dict[mode]:
        expected_sequences = list(expected_index_dict[mode][index].keys())
        for c in itertools.combinations(list(range(len(expected_sequences[0]))), hamming_distance):
            for sequence in expected_sequences:
                alternative_sequence = list(sequence)
                for sequence_index in c:
                    alternative_sequence[sequence_index] = '[ATGCN]'
                alternative_sequence = "".join(alternative_sequence)
                for instance_alternative_sequence in exrex.generate(alternative_sequence):
                    expected_index_dict[mode][index][instance_alternative_sequence] = sequence
    
    return expected_index_dict


def generate_output_file_name(output_folder:str, experiment_name:str, mode:str, index_read_num:str, fail:bool=False):
    """
    Generates output file name

    Parameters:
    -----------
    output_folder : str
        Path to user-defined output folder.
    experiment_name : str
        User-defined experiment name.
    mode : str
        User-defined mode.
    fail : bool, default False
        Indicates whether file will be used for passing or faliing reads. True is failing reads.
    index_read_num : str
        The index/read letter and number (ie I1, I2, R1, R2).

    Returns:
    --------
    file_name : str
        Full path to new file.
    """
    # define experiment specific output folder
    # TODO: maybe move this to dictionary creation function or to main()
    experiment_output_folder = os.path.join(os.path.abspath(output_folder), experiment_name)
    if not os.path.exists(experiment_output_folder):
        os.mkdir(experiment_output_folder)
    # generate file name
    if not fail:
        # TODO: check out how experiment name will be passed here - added basename piece in case it is a full path
        file_name = os.path.join(experiment_output_folder, ".".join([experiment_name, mode, index_read_num, "fq"]))
    else:
        file_name = os.path.join(experiment_output_folder, ".".join([experiment_name, "fail", index_read_num, "fq"]))
    return file_name


# TODO: make sure naming conventions match adey unidex
def generate_output_file_dict(mode_dict:dict, experiment_name:str, output_folder:str) -> dict:
    """
    Generates and opens output file locations and stores objects in dictionary

    Parameters:
    -----------
    mode_dict : dict
        Mode dictionary containing modes as keys and mode characateristics as values.
    experiment_name : str
        User-defined name of the experiment.
    output_folder : str
        Path to user-defined output folder.

    Returns:
    --------
    output_file_dict : dict
        Dictionary containing open output file objects.
    """
    # instantiate new dict
    output_file_dict = {}

    # loop through each mode and add to new dict
    for mode in mode_dict:
        output_file_dict[mode] = {
            # TODO: address single-end sequencing situations - need to be more dynamic here
            'R1_pass': open(generate_output_file_name(
                output_folder = output_folder,
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'R1'), "w"),
            'R2_pass': open(generate_output_file_name(
                output_folder = output_folder,
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'R2'), "w"),
            'R1_fail': open(generate_output_file_name(
                output_folder = output_folder,
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'R1',
                fail = True), "w"),
            'R2_fail': open(generate_output_file_name(
                output_folder = output_folder,
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'R2',
                fail = True), "w"),
            'I1_fail': open(generate_output_file_name(
                output_folder = output_folder, 
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'I1',
                fail = True), "w"),
            'I2_fail': open(generate_output_file_name(
                output_folder = output_folder,
                experiment_name = os.path.dirname(experiment_name),
                mode = mode,
                index_read_num = 'I2',
                fail = True), "w")
        }
    return output_file_dict


def close_all_files(output_file_dict:dict) -> None:
    """
    Closes all files opened in the output file dictionary.

    Parameters:
    -----------
    output_file_dict : dict
        Output file dict with keys and modes and values as subdictionaries with keys as read and index pass and value as open file.
    """
    for mode in output_file_dict:
        for instance in output_file_dict[mode]:
            output_file_dict[mode][instance].close()
    return None


def consume_next_read(open_input_file:str) -> list:
    """
    Brings next read into program.

    Parameters:
    -----------
    open_input_file : type
        Open file instance.

    Returns:
    --------
    read_fastq_components : list
        Four fastq lines for the next read.
    """
    read_fastq_components = [open_input_file.readline() for _ in range(4)]
    return read_fastq_components


def prepend_barcode_to_qname(read:list, index1:str, index2:str, index3:str) -> list:
    """
    Appends barcode to reads containing valid indexes.

    Parameters:
    -----------
    read : list
        List containing four components of read from fastq file.
    index1 : str
        The expected i7 index (already corrected if necessary).
    index2 : str
        The expected i5 index (already corrected if necessary).
    index3 : str
        The expected Tn5 index (already corrected if necessary).

    Returns:
    --------
    read : list
        Returns the same list as input, but with barcode appended between '@' and rest of sequence.
        
    Example:
    --------
    before: ['@SOMETHING:SOMETHING\n', 'AAA\n', '+\n', 'JJJ\n']
    after: ['@BARCODE:SOMETHING:SOMETHING\n', 'AAA\n', '+\n', 'JJJ\n']
    """
    # TODO: see if we need to address rev comp i7 in barcode creation
    # generate barcode
    barcode:str = "".join([index1, index2, index3])
    # prepend to qname
    read[0] = "".join(["@", barcode, ":", read[0][1:]])
    return read


def slice_read(read:list, slice:int) -> list:
    """
    Slices read to remove unwanted indexes or sequences.

    Parameters:
    -----------
    read : list
        List containing four components of read from fastq file.
    slice : int
        Number of nucleotides to remove from beginning of read.

    Returns:
    --------
    read : list
        List containing four components of read after slicing.
    """
    read[1] = read[1][slice:] # slice read sequence
    read[3] = read[3][slice:] # slice quality scores
    return read


def parse_fastq_input(
    input_files:tuple,
    mode_dict:dict,
    expected_index_dict:dict,
    output_file_dict:dict,
) -> tuple:
    """
    Processes fastqs, separating reads by mode, pass, and fail.

    Parameters:
    -----------
    input_files : tuple
        Tuple of input read and index files
    mode_dict : dict
        Mode dictionary containing modes as keys and mode characateristics as values.
    expected_index_dict : dict
        Dictionary containing modes as keys with values as subdictionaries containing index name
        as key and a list of expected indexes as values.
    output_file_dict : dict
        Dictionary containing open output file objects.

    Returns:
    --------
    total_reads : int
        Total reads processed.
    passed_reads : int
        Total reads that passed and were output to passing files.
    failed_reads : int
        Total reads that failed and were output to failing files.
    corrected_barcodes: int
        Total barcodes that were succesfully corrected.
    """
    # instantiate counters
    total_reads:int = 0
    passed_reads:int = 0
    failed_reads:int = 0
    corrected_barcodes:int = 0

    # open all input files
    open_input_files:list = [gzip.open(input_file, "rt") for input_file in input_files]
    
    # process input files
    while True:
        # TODO: account for single-end instances
        # consume next read
        read1_read, read2_read, index1_read, index2_read = [consume_next_read(open_input_file) for open_input_file in open_input_files]
        if not read1_read[0]: # break at end of files
            break
        total_reads += 1 # increment total reads

        mode_count:int = 0 # isntantiate mode count tracker - used to determine if all modes checked and should write reads to fail
        for mode in mode_dict:
            # TODO: make more dynamic (ie index4)
            # extract each index len and seq
            index1_len:int = mode_dict[mode]['index1']['index1']
            index1_seq:str = index1_read[1].strip()[:index1_len]
            index2_len:int = mode_dict[mode]['index2']['index2']
            index2_seq:str = index2_read[1].strip()[:index2_len]
            index3_len:int = mode_dict[mode]['read2']['index3'] # TODO: may need to make storage of index3 more dynamice
            index3_seq:str = read2_read[1][:index3_len]

            # TODO: make more dynamic (ie index4) - can just create a function that returns corrected sequence or None for each index in mode
            # check for matching indexes
            true_index1_seq:str = expected_index_dict[mode]['index1'][index1_seq] if index1_seq in expected_index_dict[mode]['index1'] else None
            true_index2_seq:str = expected_index_dict[mode]['index2'][index2_seq] if index2_seq in expected_index_dict[mode]['index2'] else None
            true_index3_seq:str = expected_index_dict[mode]['index3'][index3_seq] if index3_seq in expected_index_dict[mode]['index3'] else None

            # write to output if indexes match or are within hamming distance
            if true_index1_seq is not None and true_index2_seq is not None and true_index3_seq is not None:
                # increment corrected barcodes if at least one index was corrected
                # TODO: tracking the number of corrected barcodes is non-functional with new hash table index metho
                if true_index1_seq != index1_seq or true_index2_seq != index2_seq or true_index3_seq != index3_seq:
                    corrected_barcodes += 1 # increment corrected barcodes since all the way through corrections
                # prepend barcode to qname of each read
                read1_read, read2_read = [prepend_barcode_to_qname(read, true_index1_seq, true_index2_seq, true_index3_seq) for read in (read1_read, read2_read)]
                # slice reads if necessary
                if sum([mode_dict[mode]['read2'][key] for key in mode_dict[mode]['read2']]) > 0:
                    read2_read = slice_read(read2_read, sum([mode_dict[mode]['read2'][key] for key in mode_dict[mode]['read2']]))
                # write reads to passing output files
                output_file_dict[mode]['R1_pass'].write("".join(read1_read))
                output_file_dict[mode]['R2_pass'].write("".join(read2_read))
                
                passed_reads += 1 # count the passed read
                break # break out of mode_dict loop

            mode_count += 1 # increment mode count

            # if all modes checked and no pass then write to fail
            if mode_count == len(mode_dict):
                output_file_dict[mode]['R1_fail'].write("".join(read1_read))
                output_file_dict[mode]['R2_fail'].write("".join(read2_read))
                output_file_dict[mode]['I1_fail'].write("".join(index1_read))
                output_file_dict[mode]['I2_fail'].write("".join(index2_read))
                failed_reads += 1 # count the failed read
        
        # update statment
        if total_reads % 1000000 == 0:
            logging.info("{} reads processed".format(total_reads))

    # close input files
    for open_input_file in open_input_files:
        open_input_file.close()

    return total_reads, passed_reads, failed_reads, corrected_barcodes


def output_summary(
    total_reads:int,
    passed_reads:int,
    failed_reads:int,
    corrected_barcodes:int
) -> None:
    """
    Outputs summary stats to console.

    Parameters:
    -----------
    total_reads : int
        Total reads processed.
    passed_reads : int
        Total reads that passed and were output to passing files.
    failed_reads : int
        Total reads that failed and were output to failing files.
    corrected_barcodes: int
        Total barcodes that were succesfully corrected.
    """
    logging.info("Total reads processed: {}".format(total_reads))
    logging.info("Total passing reads: {}".format(passed_reads))
    logging.info("Total failed reads: {}".format(failed_reads))
    logging.info("Total corrected barcodes: {}".format(corrected_barcodes))
    return None


# define program
def main():
    args = parse_args()

    # first check if mode file is queried
    if args.query_mode_file:
        print_available_modes(args.mode_config_file)
    
    # see if mode info is requested
    if args.request_mode_info is not None:
        print_mode_details(args.mode_config_file, args.request_mode_info)

    # validate user-specified run options
    validate_run_options(args)

    # modes processing
    if args.mode_list is None:
        sys.exit("ERROR: Modes list must be specified with one or more modes!")
    mode_list:list = parse_comma_separated_inputs(
        comma_separated_input_string = args.mode_list
    )
    mode_dict:dict = generate_mode_dict(mode_list, args.mode_config_file)

    # TODO: still need to address how annotation file fits into all of this
    # annot file processing
    if args.annotation_files is not None:
        annotation_files_list:list = parse_comma_separated_inputs(
            comma_separated_input_string = args.annotation_files
        )
        annotation_dict:dict = generate_annotation_dict(
            annotation_files_list = annotation_files_list,
            mode_list = mode_list
        )
    
    # TODO: build this out
    # delayed mode functionality
    if args.delayed_mode:
        execute_delayed_mode()

    # TODO: need to address single-end instances
    # define read and index input files
    input_files:tuple = define_input_files(args)

    # generate expected index dictionary
    expected_index_dict:dict = generate_expected_index_dict(
        mode_dict = mode_dict,
        hamming_distance = args.max_hamming_distance
    )

    # generate and open output file objects
    output_file_dict:dict = generate_output_file_dict(
        mode_dict = mode_dict,
        experiment_name = args.run_folder,
        output_folder = args.output_folder
    )
    
    # parse fastq input
    total_reads, passed_reads, failed_reads, corrected_barcodes = parse_fastq_input(
        input_files = input_files,
        mode_dict = mode_dict,
        expected_index_dict = expected_index_dict,
        output_file_dict = output_file_dict
    )

    # close all files
    close_all_files(output_file_dict)

    # log results
    output_summary(total_reads, passed_reads, failed_reads, corrected_barcodes)


# run prgram
if __name__ == "__main__":
    logging.basicConfig(
        format = '%(asctime)s: %(levelname)s: %(message)s',
        level = logging.INFO
    )
    main()
