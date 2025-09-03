#!/usr/bin/env python3
"""
Convert CRYOFF .off force field files to OpenMM .xml format.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Convert CRYOFF .off force field files to OpenMM .xml format"
    )
    parser.add_argument("input_file", help="Input .off file path")
    parser.add_argument("output_file", help="Output .xml file path")
    parser.add_argument(
        "-molnames", 
        help="Specify which molecules to include in the .xml force field file"
    )
    parser.add_argument(
        "-charges", 
        help="File with Atom Name and Charge columns for charge assignment"
    )
    
    args = parser.parse_args()
    
    # TODO: Implement force field conversion logic
    print(f"Converting {args.input_file} to {args.output_file}")
    if args.molnames:
        print(f"Including molecules: {args.molnames}")
    if args.charges:
        print(f"Using charges from: {args.charges}")


if __name__ == "__main__":
    main()