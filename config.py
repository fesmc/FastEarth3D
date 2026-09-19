"""
Generate the top-level Makefile for a given machine/compiler by inserting a
compiler-configuration fragment into the template at config/Makefile.

Usage:
    python config.py config/<machine>_<compiler>      # e.g. config/macbook_gfortran

This mirrors the "configme" approach used by CLIMBER-X and Yelmo: the template
(config/Makefile) carries the build logic and a single <COMPILER_CONFIGURATION>
placeholder; the machine fragment carries everything machine-specific (compiler,
flags, netCDF paths). The shared dependency wiring lives in config/common.mk and
the source/rule lists in config/Makefile_fastearth.mk, both included by the
template.
"""

import argparse
import os
import sys

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    "config",
    metavar="CONFIG",
    type=str,
    help="path to the compiler config fragment (e.g. config/macbook_gfortran)",
)
args = parser.parse_args()

target_dir = "./"
config_dir = "config/"
config_path = args.config

if not os.path.isfile(config_path):
    sys.exit(f"error: config fragment not found: {config_path}")

template = open(config_dir + "Makefile").read()
fragment = open(config_path).read()

if "<COMPILER_CONFIGURATION>" not in template:
    sys.exit("error: config/Makefile has no <COMPILER_CONFIGURATION> placeholder")

makefile = template.replace("<COMPILER_CONFIGURATION>", fragment)
open(target_dir + "Makefile", "w").write(makefile)

print(f"\nMakefile configuration complete for: {config_path}\n")
print("Next:")
print("    make fastearth-static   # build the libfastearth.a static library")
print("    make check              # build and run the test suite")
print("    make clean              # remove objects and binaries\n")
