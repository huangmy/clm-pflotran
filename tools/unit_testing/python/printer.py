"""Utility code for printing text.

Unfortunately, the current interface uses "print", which is a keyword in
Python 2. For now, that means that users of the "print" methods must import
the Python 3 function from __future__ in order to run with Python 2.

Public classes:
Printer - Print output with some basic formatting.
CMakePrinter - Print CMake source file, with some formatting for comments.
"""
# Python 3 compatible printing in Python 2.
from __future__ import print_function

import sys

__all__ = ("Printer", "CMakePrinter")

# Declaring this to be a descendant of "object" makes it a "new-style"
# class. Done to prevent a clash with Python 3, which makes all classes
# new-style whether you explicitly inherit from "object" or not.
class Printer(object):

    """Catch-all class for formatting output.

    Public routines:
    __init__ - Constructor.
    print - Print input arguments.
    comment - Print input arguments, as a comment.
    draw_rule - Print a separator.
    print_header - Print a section header.
    print_error - Print a string as an error message.
    """

    def __init__(self, output=sys.stdout, error=sys.stderr, color=True):
        """Initialize a Printer.

        Arguments:
        output - File object to print to.
        error - File object to print errors to.
        color - Turns colors produced using ANSI control sequences on/off.
        """
        self._output = output
        self._error = error
        self._color = color

    def print(self, item, end="\n"):
        """Print the argument.

        This is intended to be used like the builtin "print" function.
        However, for compatibility with Python 2.6, we do not use a
        variadic argument, so only one item may be printed per statement.

        Arguments:
        item - Object to be printed.
        end - String appended to the end.
        """
        self._output.write(str(item)+end)

    def comment(self, string):
        """Print the input str as a comment.

        If not printing source code, printer.comment(string) is the same as
        printer.print(string).
        """
        self.print(string)

    # Utility functions for printing to the terminal.
    def draw_rule(self, char="=", length=50):
        """Draw a horizontal line that acts as a separator.

        Arguments:
        char - Character that the line is composed of.
        length - Horizontal line length.
        """
        self.comment(char*length)

    def print_header(self, string):
        """Write a string into a header, denoting a new output section."""
        self.draw_rule()
        self.comment(string)
        self.draw_rule()

    def print_error(self, error_message):
        """Print an error message to the terminal."""
        import curses.ascii

        if self._color:
            # ANSII sequence turns the text bright red.
            esc_char = chr(curses.ascii.ESC)
            to_red_text = esc_char+"[1;31m"
            to_default_text = esc_char+"[0m"
        else:
            to_red_text = ""
            to_default_text = ""

        self._error.write(to_red_text+"ERROR: "+
                          error_message+to_default_text+"\n")


class CMakePrinter(Printer):

    """Specialized printer object for writing a CMake Macros file.

    Modified public routines:
    __init__ - Restricts constructor (no longer takes a "color" argument).
    comment - Overrides parent routine: print a CMake comment.
    print_header - Extends parent routine: prints a section header.
    """

    def __init__(self, output=sys.stdout, error=sys.stderr):
        """Initialize a CMakePrinter.

        Arguments:
        output - File object to print to.
        error - File object to print errors to.

        The intention is that output is the file to push CMake output to,
        but error still prints log messages (e.g. to the terminal).
        """
        super(CMakePrinter, self).__init__(output, error, color=False)

    def comment(self, item):
        """Write a CMake comment (prepends "#")."""
        self.print("# "+item)

    def print_header(self, string):
        """Write a header in a CMake comment.

        Only difference from vanilla printer is that it adds a bit more
        space for readability.
        """
        self.print("")
        super(CMakePrinter, self).print_header(string)
        self.print("")
