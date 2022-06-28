#!/usr/bin/env python

import glob
import os
import pathlib
import pickle
import sys


def main():
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
