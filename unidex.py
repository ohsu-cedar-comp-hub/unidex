#!/usr/bin/env python

"""
Unidex

Unidex is a universal demultiplexer that can be applied to a range of sequencing data types. The
progam creates demultiplexed fastq files given gzipped fastq files, expected sequences, and 
annotations as input. Unidex was originally generated in Perl by Dr. Andrew Adey at Oregon Health
and Science University (OHSU). The tool has been translated into Python by James Adler at the Cancer
Early Detection Advanced Research Center (CEDAR), a part of the Knight Cancer Research Institue
at OHSU.
"""

__author__ = "James Adler"
__copyright__ = ""
__credits__ = ["James Adler", "Andrew Adey", "Andrew Nishida"]
__license__ = ""
__version__ = ""
__maintainer__ = "James Adler"
__email__ = "adlerj@ohsu.edu"
__status__ = "Development"


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
import time
from operator import itemgetter


def parse_args():
    parser = argparse.ArgumentParser(description = "Demultiplex fastq")

    # Info options
    info = parser.add_argument_group("Info options")
    info.add_argument("-L", "--query_mode_file", action = 'store_true', help = "List modes present in the mode config file."
                            "Can specify a different mode file with -m and it will list modes in that file."
                            "Can also provide an argument to match to refine list, e.g. 's3'")
    info.add_argument("-I", "--request_mode_info", help = "Provide info on one or more comma separated modes as detailed in the specified"
                            "modes file (-m or default).")

    # Run options
    run = parser.add_argument_group("Run options")
    run.add_argument("-R", "--run_folder", help = "Run Folder (where fastq files are present)")
    run.add_argument("-M", "--mode_list", type = str, help = "Mode list - modes must be specified in the modes.cfg file."
                            "Modes must be comma separated and will demultiplex in specified order listed.")
    run.add_argument("-l", "--delayed_mode", action = 'store_true', help = "Delayed mode. Will wait until fastq files are propagated"
                            "in the specified fastq directory (-r), then will run."
                            "Only works when specifying run name, not individual fastq files.")

    # Default options
    default = parser.add_argument_group("Default options")
    default.add_argument("-O", "--output_folder", help = "Output folder (def = run name, -R)")
    default.add_argument("-d", "--max_hamming_distance", help = "Max allowed hamming distance", type = int, default = 2)

    # Default locations
    default_locations = parser.add_argument_group("Default locations")
    default_locations.add_argument("-r", "--fastq_folder", help = "Fastq folder full path.", type = str)
    default_locations.add_argument("-m", "--mode_config_file", help = "Mode config file", required = True)

    # Fastq input (default = auto detect):
    fastq_input = parser.add_argument_group("Fastq input")
    fastq_input.add_argument("-1", "--read1_file", help = "Read 1 fastq")
    fastq_input.add_argument("-4", "--read2_file", help = "Read 2 fastq")
    fastq_input.add_argument("-2", "--index1_file", help = "Index 1 fastq")
    fastq_input.add_argument("-3", "--index2_file", help = "Index 2 fastq")

    # Other options
    other = parser.add_argument_group("Other options")
    other.add_argument("-A", "--annotation_files", type = str, help = "Annotation file(s), comma separated with mode specified"
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
    if args.fastq_folder is None:
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

    logging.info("Generating mode dictionary.")

    mode_dict:dict = {} # instantiate mode dictionary
    modes_added:int = 0 # instantiate mode tracker

    # read through config file and store entries specified in modes list
    open_mode_config_file = open(mode_config_file, 'r')

    while modes_added < len(mode_list):
        line = open_mode_config_file.readline()
        if not line: # break at end of file
            break
        if line.startswith("#") or not line[0]: # skip header / description lines
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

    
    ################ DEPRECATED BUT SAVING FOR NOW ##################
    # # extract components of mode
    # mode, read1, index1, index2, read2, index_files = mode_components[0], mode_components[1], mode_components[2], mode_components[3], mode_components[4], mode_components[5:]
    ################ DEPRECATED BUT SAVING FOR NOW ##################



    # instantiate mode and compoenents within mode dictionary
    mode = mode_components.pop(0)
    mode_dict[mode] = {} # add mode
    while '=' not in mode_components[0]:
        designation = mode_components.pop(0)
        split_designation = designation.split(',')
        components_encountered_in_designation:list = []
        # track start and end position
        start_pos:int = 0
        for component in split_designation:
            component_id, component_len = component.split(':')
            component_len = int(component_len)
            if component_id in mode_dict[mode]:
                sys.exit("\nERROR: {} is specified in multiple locations in mode file!!!\n")
            if component_id != 'null':
                print(designation)
                mode_dict[mode][component_id] = {
                    'is_index': True if 'index' in component_id else False,
                    'index_in_read': True if 'index' in designation and 'read' in designation else False,
                    'start_pos': start_pos,
                    'end_pos': start_pos + component_len,
                    'needs_trimmed': True if 'null' in designation and 'read' in component_id or 'index' in designation and 'read' in component_id else False,
                    'trim_len': 0,
                    'location': 'read1' if 'read1' in designation else 'read2' if 'read2' in designation else component_id # this will be 'index1', 'index2', 'read1', or 'read2' to be referenced in parse
                }
                components_encountered_in_designation.append(component_id)
            # increment start position
            start_pos += component_len
        # have to go back at end and update trim length for reads - which is euqual to the up to date start position
        # TODO: this may not be the case if there are also indexes at the end of the read or in middle of read
        # TODO: could loop back through the components and extract the lengths to get this right in the future
        for component_id in components_encountered_in_designation:
            if 'read' in component_id:
                mode_dict[mode][component_id]['trim_len'] = start_pos



    ########## DEPRECATED BUT SAVING FOR NOW ##################        
    # mode_dict[mode] = {}
    # mode_dict[mode]['read1'] = {part.split(":")[0]: int(part.split(":")[1]) for part in read1.split(",")}
    # mode_dict[mode]['index1'] = {part.split(":")[0]: int(part.split(":")[1]) for part in index1.split(",")}
    # mode_dict[mode]['index2'] = {part.split(":")[0]: int(part.split(":")[1]) for part in index2.split(",")}
    # mode_dict[mode]['read2'] = {part.split(":")[0]: int(part.split(":")[1]) for part in read2.split(",")}
    ########## DEPRECATED BUT SAVING FOR NOW ################## 



    # add index file paths to dictionary - ignoring special flags
    mode_dict[mode]['index_file_paths'] = {}
    while mode_components:
        index_path_component:str = mode_components.pop(0)
        index_id, index_path = index_path_component.split('=')
        mode_dict[mode]['index_file_paths'][index_id] = index_path
    print(mode_dict)
    
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
        Dictionary to be used for mapping barcodes to annotation subjects.
    annotation_subjects_dict : dict
        Dictionary to be used for generating output files.
    """
    annotation_dict:dict = {} # instantiate empty annotation dictionary
    annotation_subjects_dict:dict = {} # instantiate empty dictionary of unique subjects by mode

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
        if mode not in annotation_subjects_dict:
            annotation_subjects_dict[mode] = set()
        annotation_dict[mode][cellID] = annot
        annotation_subjects_dict[mode].add(annot)

        # process entire file
        while True:
            line = open_annotation_file.readline().strip().split()            
            if not line: # break if end of file
                break
            if mode_in_annot_file: # may be multiple modes within annotation file
                mode, cellID, annot = line # extract all three columns
                if mode not in annotation_dict:
                    annotation_dict[mode] = {}
                if mode not in annotation_subjects_dict:
                    annotation_subjects_dict[mode] = set()
            else:
                cellID, annot = line # extract only two columns if mode specified in command line        
            annotation_dict[mode][cellID] = annot # add instance to dictionary
            annotation_subjects_dict[mode].add(annot) # add instance to subject dictionary
            
        logging.info("Total lines processed for annotation file '{}': {}".format(annotation_file, len(annotation_dict[mode])))
        open_annotation_file.close()

    return annotation_dict, annotation_subjects_dict


# TODO: build this out
def execute_delayed_mode():
    sys.exit("Delayed mode not currently implemented.")


def define_input_files(args) -> tuple[str, str, str, str]:
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
    elif args.fastq_folder is not None:
        read1_file = os.path.abspath(os.path.join(args.fastq_folder, "Undetermined_S0_L001_R1_001.fastq.gz"))
        if not os.path.exists(read1_file):
            sys.exit("ERROR: Read1 file does not exist {}".format(read1_file))
    else:
        sys.exit("ERROR: Either read1 file (-1 flag) OR run folder (-R flag) must be defined.")

    # define read 2
    if args.read2_file is not None:
        read2_file = os.path.abspath(args.read2_file)
    elif args.fastq_folder is not None:
        read2_file = os.path.abspath(os.path.join(args.fastq_folder, "Undetermined_S0_L001_R2_001.fastq.gz"))
        if not os.path.exists(read2_file):
            logging.info("No read 2 file detected at path {}.\nMoving forward with single end read.".format(read2_file))
            read2_file = ""

    # define index 1
    if args.index1_file is not None:
        index1_file = os.path.abspath(args.index1_file)
    elif args.fastq_folder is not None:
        index1_file = os.path.abspath(os.path.join(args.fastq_folder, "Undetermined_S0_L001_I1_001.fastq.gz"))
        if not os.path.exists(index1_file):
            logging.info("No index 1 file detected at path {}".format(index1_file))
            index1_file = ""

    # define index 2    
    if args.index2_file is not None:
        index2_file = os.path.abspath(args.index2_file)
    elif args.fastq_folder is not None:
        index2_file = os.path.abspath(os.path.join(args.fastq_folder, "Undetermined_S0_L001_I2_001.fastq.gz"))
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
    hamming_distance : 2
        Maximum allowable hamming distance for each index.

    Returns:
    --------
    expected_index_dict : dict
        Dictionary containing modes as keys with values as subditionaries containing hamming distance
        possible indexes as keys and true indexes as values.
    """

    logging.info("Generating dictionary of expected indexes.")

    # instantiate dicts
    expected_index_dict:dict = {} # dictionary mapping hamming distance indexes to true indexes

    # extract modes and indexes from mode_dict and defined index files
    for mode in mode_dict:
        expected_index_dict[mode] = {}
        for index in mode_dict[mode]['index_file_paths']:
            with open(mode_dict[mode]['index_file_paths'][index], "r") as f:
                expected_index_dict[mode][index] = {line.split('\t')[-1]:line.split('\t')[-1] for line in f.read().split('\n') if line}
    
        # generate hamming distance possibilities
        ignore_hamming_distance_priority:bool = False # avoids overwriting higher priority (lower hamming distance) indexes
        if hamming_distance > 0:
            for index in expected_index_dict[mode]:
                # generate list of true indexes for 
                expected_sequences = list(expected_index_dict[mode][index].keys())
                encountered_sequences = set(expected_sequences)
                # loop through hamming distance step wise
                for instance_hamming_distance in range(1, hamming_distance + 1):
                    logging.info("Generating alternative {} sequences for hamming distance: {}".format(index, instance_hamming_distance))
                    hamming_distance_alternative_sequences_encountered:set = set() # tracks indexes of higher priority
                    # generate index combinations for given hamming distance
                    for c in itertools.combinations(list(range(len(expected_sequences[0]))), instance_hamming_distance):
                        # convert index combinations to sequence combinations
                        for sequence in expected_sequences:
                            alternative_sequence = list(sequence)
                            for sequence_index in c:
                                alternative_sequence[sequence_index] = '[ATGCN]'
                            alternative_sequence = "".join(alternative_sequence)
                            alternative_sequences = exrex.generate(alternative_sequence)
                            # add appropriate alternative sequences to expected sequences dict
                            for instance_alternative_sequence in alternative_sequences:
                                if not instance_alternative_sequence in expected_index_dict[mode][index]:
                                    expected_index_dict[mode][index][instance_alternative_sequence] = sequence
                                    hamming_distance_alternative_sequences_encountered.add(instance_alternative_sequence) # tracks current hamming distance instances
                                elif instance_alternative_sequence not in expected_sequences and ignore_hamming_distance_priority: # avoids marking true index as ambiguous
                                    expected_index_dict[mode][index][instance_alternative_sequence] = 'ambiguous'
                                elif instance_alternative_sequence in hamming_distance_alternative_sequences_encountered and not ignore_hamming_distance_priority: # avoids removing indexes of higher priority
                                    expected_index_dict[mode][index][instance_alternative_sequence] = 'ambiguous'
                        # add hamming distance indexes to indexes list of higher priority
                        encountered_sequences.update(hamming_distance_alternative_sequences_encountered)
    return expected_index_dict


def generate_output_file_name(
    output_folder:str,
    experiment_name:str,
    mode:str,
    index_read_num:str,
    fail:bool=False,
    annotation_subject:str=None
):
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
    index_read_num : str
        The index/read letter and number (ie I1, I2, R1, R2).
    fail : bool, default False
        Indicates whether file will be used for passing or faliing reads. True is failing reads.
    annotation_subject : str, default None
        If an annotation file was specified, then output file names will need to include annotation subject name.

    Returns:
    --------
    file_name : str
        Full path to new file.
    """
    # define experiment specific output folder
    # TODO: maybe move this to dictionary creation function or to main()
    experiment_output_folder = os.path.join(os.path.abspath(output_folder), experiment_name)
    logging.info("Outputting files to: {}".format(experiment_output_folder))
    if not os.path.exists(experiment_output_folder):
        os.mkdir(experiment_output_folder)
    # generate file name
    if not fail:
        if annotation_subject is None:
            file_name = os.path.join(experiment_output_folder, ".".join([experiment_name, mode, index_read_num, "fq"]))
        else:
            file_name = os.path.join(experiment_output_folder, ".".join([experiment_name, mode, annotation_subject, index_read_num, "fq"]))
    else:
        file_name = os.path.join(experiment_output_folder, ".".join([experiment_name, "fail", index_read_num, "fq"]))
    return file_name


# TODO: make sure naming conventions match adey unidex
def generate_output_file_dict(
    mode_dict:dict,
    experiment_name:str,
    output_folder:str,
    annotation_subjects_dict:dict=None
) -> dict:
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
    annotation_subjects_dict : dict, default None
        Dictionary to be used for generating output files.

    Returns:
    --------
    passing_output_file_dict : dict
        Dictionary containing open output file objects for passing reads.
    """
    # instantiate output file dict with fail keys (not mode dependent)
    failing_output_file_dict:dict = {
        'R1_fail': open(generate_output_file_name(
            output_folder = output_folder,
            experiment_name = experiment_name,
            mode = None,
            index_read_num = 'R1',
            fail = True), "w"),
        'R2_fail': open(generate_output_file_name(
            output_folder = output_folder,
            experiment_name = experiment_name,
            mode = None,
            index_read_num = 'R2',
            fail = True), "w"),
        'I1_fail': open(generate_output_file_name(
            output_folder = output_folder, 
            experiment_name = experiment_name,
            mode = None,
            index_read_num = 'I1',
            fail = True), "w"),
        'I2_fail': open(generate_output_file_name(
            output_folder = output_folder,
            experiment_name = experiment_name,
            mode = None,
            index_read_num = 'I2',
            fail = True), "w")
    }

    # instantiate dictionary for passing output file
    passing_output_file_dict:dict = {}

    # if no annotations specified then generate simple dictionary
    if annotation_subjects_dict is None:
        # loop through each mode and add to new dict
        for mode in mode_dict:
            passing_output_file_dict[mode] = {
                # TODO: address single-end sequencing situations - need to be more dynamic here
                'R1_pass': open(generate_output_file_name(
                    output_folder = output_folder,
                    experiment_name = experiment_name,
                    mode = mode,
                    index_read_num = 'R1'), "w"),
                'R2_pass': open(generate_output_file_name(
                    output_folder = output_folder,
                    experiment_name = experiment_name,
                    mode = mode,
                    index_read_num = 'R2'), "w")
            }
    # if there are annotation files specified, the output files must include annotation names
    else:
        for mode in annotation_subjects_dict:
            passing_output_file_dict[mode] = {}
            for annotation_subject in annotation_subjects_dict[mode]:
                passing_output_file_dict[mode][annotation_subject] = {
                    'R1_pass': open(generate_output_file_name(
                        output_folder = output_folder,
                        experiment_name = experiment_name,
                        mode = mode,
                        index_read_num = 'R1',
                        annotation_subject = annotation_subject), 'w'),
                    'R2_pass': open(generate_output_file_name(
                        output_folder = output_folder,
                        experiment_name = experiment_name,
                        mode = mode,
                        index_read_num = 'R2',
                        annotation_subject = annotation_subject), 'w')
                }
    return passing_output_file_dict, failing_output_file_dict


def close_all_files(passing_output_file_dict:dict, failing_output_file_dict:dict, annotation_file_used:bool=False) -> None:
    """
    Closes all files opened in the output file dictionary.

    Parameters:
    -----------
    passing_output_file_dict : dict

    failing_output_file_dict : dict

    annotation_file_used : bool, default False
        Indicates if annotation file was used. Changes structure of the passing_output_file_dict. 

    """
    for instance in failing_output_file_dict:
        failing_output_file_dict[instance].close()
    for mode in passing_output_file_dict:
        if annotation_file_used:
            for annotation_subject in passing_output_file_dict[mode]:
                for instance in passing_output_file_dict[mode][annotation_subject]:
                    passing_output_file_dict[mode][annotation_subject][instance].close()
        else:
            for instance in passing_output_file_dict[mode]:
                passing_output_file_dict[mode][instance].close()
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


def prepend_barcode_to_qname(read:list, true_index_seqs:list) -> list:
    """
    Appends barcode to reads containing valid indexes.

    Parameters:
    -----------
    read : list
        List containing four components of read from fastq file.
    true_index_seqs : list
        List containing all true indexes.

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
    barcode:str = "".join(true_index_seqs)
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
    passing_output_file_dict:dict,
    failing_output_file_dict:dict,
    annotation_dict:dict=None
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
    passing_output_file_dict : dict

    failing_output_file_dict : dict

    annotation_dict : dict, default None
        Annotation dictionary containing associations between barcodes and subjects.

    Returns:
    --------
    total_reads : int
        Total reads processed.
    passed_reads : int
        Total reads that passed and were output to passing files.
    failed_reads : int
        Total reads that failed and were output to failing files.
    corrected_barcodes : int
        Total barcodes that were succesfully corrected.
    unspecified_barcodes : int
        Total barcodes thrown out due to unspecified annotation.
    """
    # instantiate counters
    total_reads:int = 0
    passed_reads:int = 0
    failed_reads:int = 0
    corrected_barcodes:int = 0
    ambiguous_barcodes:int = 0 # tracks reads thrown out due to hamming distance collision
    unspecified_barcodes:int = 0 # tracks reads thrown out due to unspecified annotation

    # open all input files
    open_input_files:dict = {}
    for input_file in input_files:
        if '_R1_' in input_file:
            open_input_files['read1'] = gzip.open(input_file, "rt")
        elif '_R2_' in input_file:
            open_input_files['read2'] = gzip.open(input_file, "rt")
        elif '_I1_' in input_file:
            open_input_files['index1'] = gzip.open(input_file, "rt")
        elif '_I2_' in input_file:
            open_input_files['index2'] = gzip.open(input_file, "rt")
        else:
            continue
        logging.info("Opening input file: {}".format(input_file))
        
    # process input files
    while True:
        # TODO: account for single-end instances
        # consume next read
        reads:dict = {}
        for read in open_input_files:
            reads[read] = consume_next_read(open_input_files[read])
        if not reads['read1'][0]: # break at end of files
            break
        total_reads += 1 # increment total reads

        mode_count:int = 0 # isntantiate mode count tracker - used to determine if all modes checked and should write reads to fail
        ambiguous_index_encountered:bool = False # tracks if ambiguous read found due to hamming distance collision
        unspecified_annotation:bool = False # tracks if unspecified annotation encountered
        
        for mode in mode_dict:

            mode_count += 1 # increment mode count

            # TODO: make more dynamic (ie index4) - can just make this a dictionary instead of individual objects
            # extract each index len and seq

            observed_index_seqs:list = []
            true_index_seqs:list = []
            
            for designation in mode_dict[mode]:
                if designation != 'index_file_paths': # TODO: come back and adjust where this is stored
                    if 'index' in designation:
                        # identify location of index sequence
                        location = mode_dict[mode][designation]['location']
                        start_pos = mode_dict[mode][designation]['start_pos']
                        end_pos = mode_dict[mode][designation]['end_pos']
                        # extract index sequence of interest from read
                        read_sequence = reads[location][1]
                        index_sequence = read_sequence[start_pos:end_pos]
                        observed_index_seqs.append(index_sequence)
                        # determine true index sequence
                        if index_sequence in expected_index_dict[mode][designation]:
                            true_index_seqs.append(expected_index_dict[mode][designation][index_sequence])
                        else:
                            true_index_seqs.append(None)



            ################ DEPRECATED BUT SAVING FOR NOW ########################                        
            # index1_len:int = mode_dict[mode]['index1']['index1']
            # index1_seq:str = reads['index1'][1].strip()[:index1_len]
            # index2_len:int = mode_dict[mode]['index2']['index2']
            # index2_seq:str = reads['index2'][1].strip()[:index2_len]
            # index3_len:int = mode_dict[mode]['read2']['index3'] # TODO: may need to make storage of index3 more dynamic
            # index3_seq:str = reads['read2'][1][:index3_len]                     
            
            # # TODO: make more dynamic (ie index4) - can just create a function that returns corrected sequence or None for each index in mode
            # # check for matching indexes
            # true_index1_seq:str = expected_index_dict[mode]['index1'][index1_seq] if index1_seq in expected_index_dict[mode]['index1'] else None
            # true_index2_seq:str = expected_index_dict[mode]['index2'][index2_seq] if index2_seq in expected_index_dict[mode]['index2'] else None
            # true_index3_seq:str = expected_index_dict[mode]['index3'][index3_seq] if index3_seq in expected_index_dict[mode]['index3'] else None
            ################ DEPRECATED BUT SAVING FOR NOW ########################   



            # if all legitimate indexes
            if None not in true_index_seqs:
                # hamming distance collision
                if 'ambiguous' in true_index_seqs:
                    ambiguous_index_encountered = True
                # write to output if indexes match or are within hamming distance and no hamming distance collision
                else:
                    # increment corrected barcodes if at least one index was corrected
                    if observed_index_seqs != true_index_seqs:
                        corrected_barcodes += 1 # increment corrected barcodes since all the way through corrections
                    # prepend barcode to qname of each read
                    read1_read, read2_read = [prepend_barcode_to_qname(read, true_index_seqs) for read in itemgetter('read1', 'read2')(reads)]
                    # slice reads if necessary
                    for designation in mode_dict[mode]:
                        if 'read' in designation:
                            reads[designation] = slice_read(reads[designation], mode_dict[mode][designation]['trim_len'])
                    # write reads to passing output files
                    if annotation_dict is not None:
                        # TODO: this should probbly be a function - will need to be more dynamic
                        read_barcode = "".join(true_index_seqs)
                        # asses if barcode specified in annotation by user
                        if read_barcode in annotation_dict[mode]:
                            annotation_subject = annotation_dict[mode][read_barcode]
                            unspecified_annotation:bool = False # used to catch barcodes not specified in annotation                        
                        else:
                            logging.info("Expected barcode found from sequence but annotation not specificed.")
                            unspecified_annotation = True # used to catch barcodes not specified in annotation 
                        if not unspecified_annotation:
                            passing_output_file_dict[mode][annotation_subject]['R1_pass'].write("".join(reads['read1']))
                            passing_output_file_dict[mode][annotation_subject]['R2_pass'].write("".join(reads['read2']))
                    else:
                        passing_output_file_dict[mode]['R1_pass'].write("".join(reads['read1']))
                        passing_output_file_dict[mode]['R2_pass'].write("".join(reads['read2']))
                    if not unspecified_annotation:
                        passed_reads += 1 # count the passed read
                        break # break out of mode_dict loop

            # if all modes checked and no pass then write to fail
            if mode_count == len(mode_dict):
                failing_output_file_dict['R1_fail'].write("".join(reads['read1']))
                failing_output_file_dict['R2_fail'].write("".join(reads['read2']))
                failing_output_file_dict['I1_fail'].write("".join(reads['index1']))
                failing_output_file_dict['I2_fail'].write("".join(reads['index2']))
                failed_reads += 1 # count the failed read
                
                # count ambiguous barcode if hamming distance collision encountered
                if ambiguous_index_encountered:
                    ambiguous_barcodes += 1
                # counts barcodes thrown out due to unspecified annotation
                if unspecified_annotation:
                    unspecified_barcodes += 1
        
        # update statment
        if total_reads % 1000000 == 0:
            logging.info("{} reads processed".format(total_reads))

    # close input files
    for open_input_file in open_input_files:
        open_input_files[open_input_file].close()

    return total_reads, passed_reads, failed_reads, corrected_barcodes, ambiguous_barcodes, unspecified_barcodes


def output_summary(
    total_reads:int,
    passed_reads:int,
    failed_reads:int,
    corrected_barcodes:int,
    ambiguous_barcodes:int,
    unspecified_barcodes:int,
    output_folder:str,
    experiment_name:str
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
    ambiguous_barcodes : int
        Total barcodes thrown out due to hamming distance collision.
    unspecified_barcodes : int
        Total barcodes thrown out due to unspecified annotation.
    output_folder : str
        Folder to output summary file to.
    experiment_name : str
        Run name.
    """
    logging.info("Total reads processed: {}".format(total_reads))
    logging.info("Total passing reads: {}".format(passed_reads))
    logging.info("Total failed reads: {}".format(failed_reads))
    logging.info("Total corrected barcodes: {}".format(corrected_barcodes))
    logging.info("Total barcodes thrown out due to hamming distance collision: {}".format(ambiguous_barcodes))
    logging.info("Total barcodes thrown out due to unspecified annotation: {}".format(unspecified_barcodes))
    
    # define output summary file
    timestr = time.strftime("%Y%m%d-%H%M%S")
    output_filename = ".".join(["summary-output", timestr, "txt"])
    output_filepath = os.path.join(os.path.abspath(output_folder), experiment_name, output_filename)
    
    with open(output_filepath, 'w') as f:
        f.write("Total reads processed: {}\n".format(total_reads))
        f.write("Total passing reads: {}\n".format(passed_reads))
        f.write("Total failed reads: {}\n".format(failed_reads))
        f.write("Total corrected barcodes: {}\n".format(corrected_barcodes))
        f.write("Total barcodes thrown out due to hamming distance collision: {}\n".format(ambiguous_barcodes))
        f.write("Total barcodes thrown out due to unspecified annotation: {}\n".format(unspecified_barcodes))
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
        sys.exit("\nERROR: Modes list must be specified with one or more modes!\n")
    mode_list:list = parse_comma_separated_inputs(
        comma_separated_input_string = args.mode_list
    )
    mode_dict:dict = generate_mode_dict(mode_list, args.mode_config_file)

    # annot file processing
    if args.annotation_files is not None:
        annotation_files_list:list = parse_comma_separated_inputs(
            comma_separated_input_string = args.annotation_files
        )
        annotation_dict, annotation_subjects_dict = generate_annotation_dict(
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

    # define experiment name
    if args.run_folder is not None:
        experiment_name = os.path.basename(args.run_folder)
    else:
        experiment_name = os.path.basename(args.fastq_folder)
    logging.info("User has defined experiment name: {}".format(experiment_name))

    # generate and open output file objects
    passing_output_file_dict, failing_output_file_dict = generate_output_file_dict(
        mode_dict = mode_dict,
        experiment_name = experiment_name,
        output_folder = args.output_folder,
        annotation_subjects_dict = annotation_subjects_dict if args.annotation_files is not None else None
    )
    
    # parse fastq input
    total_reads, passed_reads, failed_reads, corrected_barcodes, ambiguous_barcodes, unspecified_barcodes = parse_fastq_input(
        input_files = input_files,
        mode_dict = mode_dict,
        expected_index_dict = expected_index_dict,
        passing_output_file_dict = passing_output_file_dict,
        failing_output_file_dict = failing_output_file_dict,
        annotation_dict = annotation_dict if args.annotation_files is not None else None
    )

    # close all files
    close_all_files(
        passing_output_file_dict = passing_output_file_dict,
        failing_output_file_dict = failing_output_file_dict,
        annotation_file_used = True if args.annotation_files is not None else False
    )

    # log results
    output_summary(
        total_reads,
        passed_reads,
        failed_reads,
        corrected_barcodes,
        ambiguous_barcodes,
        unspecified_barcodes,
        args.output_folder,
        experiment_name
    )


# run prgram
if __name__ == "__main__":
    logging.basicConfig(
        format = '%(asctime)s: %(levelname)s: %(message)s',
        level = logging.INFO
    )
    main()
