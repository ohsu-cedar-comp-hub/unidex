#!/usr/bin/env python

# unidex initial python build

# import modules
import argparse
import sys
import os
import logging
import json


def parse_args():
    parser = argparse.ArgumentParser(description = "Demultiplex fastq")
    parser.add_argument("-M", "--mode_list", type = str, help = "Mode list - modes must be specified in the modes.cfg file\n"
                            "Modes must be comma separated and will demultiplex in specified order listed.\n")
    parser.add_argument("-A", "--annotation_files", type = str, help = "Annotation file(s), comma separated with mode specified\n"
                            "If only one mode is specified, then it will default to that mode\n"
                            "[mode1]=[annot_file1],[mode2]=[annot_file2],etc... OR\n"
                            "First column of annot file designated mode for that annot")
    parser.add_argument("-m", "--mode_config_file", help = "Mode config file")
    return parser.parse_args()


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
                    "Please validate modes are specified in config file: {}".format(mode_list, mode_dict, os.path.realpath(mode_config_file)))


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
    mode, read1, index1, index2, read2, index_files = mode_components[0], mode_components[1], mode_components[2], mode_components[3], mode_components[4], mode_components[5:]
    mode_dict[mode] = {}
    mode_dict[mode]['read1'] = read1
    mode_dict[mode]['index1'] = index1
    mode_dict[mode]['index2'] = index2
    mode_dict[mode]['read2'] = read2
    mode_dict[mode]['index_files'] = index_files
    return mode_dict


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

    logging.info("Mode dicitionary successfully created:") # log dictionary
    for mode in mode_dict:
        logging.info("{}".format(mode))
        for component in mode_dict[mode]:
            logging.info("{}: {}".format(component, mode_dict[mode][component]))

    return mode_dict

def generate_annotation_dict(annotation_files_list:list, mode_list:list) -> dict:
    """
    Process annotation file(s).

    Parameters:
    -----------
    annotation_files_list : list
        List of all annotation files.
    mode_list : list
        List of corresponding mode designation for each annotation file.
    
    Returns:
    --------
    annotation_dict : dict
    """
    # instantiate empty annotation dictionary
    annotation_dict:dict = {}

    # process each annotation file
    for i, annotation_file in enumerate(annotation_files_list):
        logging.info("Processing annotation file: {}".format(annotation_file))
        annotation_dict[annotation_file] = {} # set mode to corresponding value from mode_list
        annotation_dict[annotation_file][mode_list[i]] = {}
        open_annotation_file = open(annotation_file, 'r')
        
        # process each line of annotation file
        while True:
            line = open_annotation_file.readline().strip().split()            
            if not line: # break if end of file
                break
            cellID, annot = line # exctract cellid and annotation from current line            
            annotation_dict[annotation_file][mode_list[i]][cellID] = annot # add instance dictionary
        
        logging.info("Total lines processed for annotation file '{}': {}".format(annotation_file, len(annotation_dict[annotation_file][mode_list[i]])))
        open_annotation_file.close()

    return annotation_dict


# define program
def main():
    args = parse_args()

    # modes processing
    mode_list:list = parse_comma_separated_inputs(
        comma_separated_input_string = args.mode_list
    )
    mode_dict:dict = generate_mode_dict(mode_list, args.mode_config_file)

    # annot file processing
    annotation_files_list:list = parse_comma_separated_inputs(
        comma_separated_input_string = args.annotation_files
    )
    annotation_dict:dict = generate_annotation_dict(
        annotation_files_list = annotation_files_list,
        mode_list = mode_list
    )
    print(annotation_dict)


# run prgram
if __name__ == "__main__":
    logging.basicConfig(
        format = '%(asctime)s: %(levelname)s: %(message)s',
        level = logging.INFO
    )
    main()
