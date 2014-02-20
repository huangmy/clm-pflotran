"""Functions to set up machine-specific settings.

Public functions:
add_compiler_env - Adds compile/machine settings to the environment.
"""

# Python 3 compatible printing in Python 2.
from __future__ import print_function

from os import environ
import platform
import re

from xml_utils import best_match, all_matches

__all__ = ("add_compiler_env")

def get_machine_name():
    """Returns a canonical version of the machine name."""
    name = platform.node()
    # Strip everything after the first ".", and whitespace.
    name = name.split(".")[0].strip()
    if re.match("^yslogin[0-9]+", name):
        name = "yellowstone"
    elif re.match("^miralac[0-9]+", name):
        name = "mira"
    return name

def load_machine_env(compiler):
    """Add machine environment variables not in config_compilers.xml

    Besides simply setting variables, this may also load some modules.
    """

    mach = get_machine_name()

    if mach == "yellowstone":
        environ["MODULESHOME"] = "/glade/apps/opt/lmod/lmod/"
        execfile(environ["MODULESHOME"]+"/init/env_modules_python.py",
                 globals())
        module("purge")
        module("load", "ncarenv/1.0")
        module("load", "ncarbinlibs/1.0")
        if compiler == "intel":
            module("load", "intel/13.1.2")
        module("load", "ncarcompilers/1.0")
        module("load", "cmake/2.8.10.2")


def make_cmake_macros(xml_tree, machine_dict, macros_file):
    """Write out a "Macros" file for CMake.

    Given an xml tree from config_compilers, a dictionary with machine
    options, and an open file, write a CMake macros file.
    """
    from itertools import chain

    from printer import CMakePrinter

    # Print header to file.
    macros_printer = CMakePrinter(macros_file)
    macros_printer.comment("CESM build flags for:")
    macros_printer.comment("  Compiler = "+machine_dict["COMPILER"]+":")
    macros_printer.comment("  Machine = "+machine_dict["MACH"]+":")
    macros_printer.comment("  OS = "+machine_dict["OS"]+":")

    # pFUnit location if it exists.
    match = best_match(xml_tree, "compiler/PFUNIT_PATH", machine_dict)
    if match is not None:
        macros_printer.print_header("pFUnit location.")
        macros_printer.print(
            "list(APPEND CMAKE_PREFIX_PATH "+match+")"
            )

    # Normal and debug dictionaries for looking things up in
    # config_compilers.
    normal_dict = machine_dict.copy()
    normal_dict["DEBUG"] = "FALSE"

    debug_dict = machine_dict.copy()
    debug_dict["DEBUG"] = "TRUE"

    def add_formatted_flags(flags_name, format):
        """Print CMake flags using macros_printer.

        Arguments:
        flags_name - Name to search for in config_compilers (e.g. FFLAGS).
        format - Function that takes a build type and flag match, and
                 returns the string to print out.
        """

        paths = ["compiler/"+flags_name, "compiler/ADD_"+flags_name]

        # This creates an iterable over elements in config_compilers that
        # match in non-debug mode.
        normal_matches = chain.from_iterable(
            all_matches(xml_tree, path, normal_dict) for path in paths
            )
        for match in normal_matches:
            macros_printer.print(format("CESM", match))

        # Now the same for debug mode.
        debug_matches = chain.from_iterable(
            all_matches(xml_tree, path, debug_dict) for path in paths
            )
        for match in debug_matches:
            macros_printer.print(format("CESM_DEBUG", match))

    # Below, we just use a bunch of lambda functions to describe how
    # the build type and a matching element (e.g. an FFLAGS entry) are
    # turned into a CMake function call.

    macros_printer.print_header("CPP definitions.")

    add_formatted_flags("CPPDEFS", lambda b, m:
                            "add_config_definitions("+b+" "+m+")")

    macros_printer.print_header("Fortran flags.")

    add_formatted_flags("FFLAGS", lambda b, m:
                            "add_flags(CMAKE_Fortran_FLAGS_"+b+" "+m+")")

    macros_printer.print_header("C flags.")

    add_formatted_flags("CFLAGS", lambda b, m:
                            "add_flags(CMAKE_C_FLAGS_"+b+" "+m+")")

    macros_printer.print_header("Linker flags.")

    add_formatted_flags("LDFLAGS", lambda b, m:
                            "add_flags(CMAKE_EXE_LINKER_FLAGS_"+b+" "+m+")")


def add_compiler_env(compiler, compiler_xml_path, macros_file_name):
    """Add config_compilers information to the environment.

    Given a compiler and the path to a config_compilers.xml file, adds
    compiler information to the environment, as well as processing the
    compiler flags into a file for inclusion by CMake.
    """
    from xml.etree.ElementTree import ElementTree

    machine_dict = {
        "COMPILER": compiler,
        "MACH": get_machine_name(),
        "OS": platform.system(),
        }

    # Parse config_compilers for information.
    compiler_xml_tree = ElementTree()
    compiler_xml_tree.parse(compiler_xml_path)

    def add_env_var(xml_path, var_name):
        match = best_match(compiler_xml_tree, xml_path,
                           machine_dict)
        assert match is not None, "Could not determine "+var_name+ \
            " from compiler/machines xml file."
        environ[var_name] = match.strip()

    add_env_var("compiler/SFC", "FC")
    add_env_var("compiler/SCC", "CC")

    load_machine_env(compiler)

    with open(macros_file_name, "w") as macros_file:
        make_cmake_macros(compiler_xml_tree, machine_dict, macros_file)
