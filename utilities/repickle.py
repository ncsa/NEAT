#!/usr/bin/env python
#
#
#   repickle.py
#   What is its purpose?
#
#   Can take a directory containing python pickle files and repickles them
#
#   Usage: python utilities/repickle.py path_to_dir
#
#
#  Python 3 ready



import sys
import os
import pathlib
import glob
import pickle


def main():
    '''
    Repickles every pickle file in the given directory
    '''
    dir_to_repickle = pathlib.Path(sys.argv[1])

    if not dir_to_repickle.is_dir():
        print("Input is not a directory.")
        sys.exit(1)

    os.chdir(dir_to_repickle)
    for file in glob.glob("*.p"):
        data = pickle.load(open(file, 'rb'), encoding="bytes")
        pickle.dump(data, open(file, "wb"))


if __name__ == "__main__":
    main()
