"""XML search utilities.

These are based on the standard library xml.etree.ElementTree module.
However, that module is not explicitly imported at this time, so it is in
principle possible to pass arguments that imitate the interface of that
module.

Exported functions:
best_match - Search for a specific element in an XML tree.
all_matches - Search for all matching elements in an XML tree.
"""

from comparable import Comparable

__all__ = ("best_match", "all_matches")

class ElementMatch(Comparable):

    """Class to aid searching for/matching an xml element.

    Public methods:
    __init__ - Create an ElementMatch.
    __iadd__ - Add in "quality" of other match.
    __eq__ - Test if quality is same as another match.
    __lt__ - Compare quality to another match.
    __nonzero__/__bool__ - Test for a valid match.
    get_text - Get the text within the element.

    Note that using __eq__ to test quality interferes with hashable
    collections, so for now, __hash__ = None.
    """

    def __init__(self, element=None, quality=0):
        """Define a element match.

        Arguments:
        element - The element that was found.
        quality - Abstract measure of "quality" (specificity) of the match.

        With both options left out, this will return a null ElementMatch.
        """
        self._element = element
        self._quality = quality

    def __iadd__(self, other):
        """Add the quality from another match to that of this match."""
        self._quality += other._quality
        return self

    __hash__ = None
    """Not implemented, because __eq__ is overridden."""

    def __eq__(self, other):
        """Check whether two matches have the same quality."""
        return (not self and not other) or \
            (self and other and self._quality == other._quality)

    def __lt__(self, other):
        """Compare matches by their quality; null matches always lose."""
        return (not self and other) or \
            (self and other and self._quality < other._quality)

    def __bool__(self):
        """Test if this is a real match or just the null one."""
        return self._element is not None

    __nonzero__ = __bool__

    def get_text(self):
        """Get text from valid matches, or None for null matches."""
        if self:
            return self._element.text

def element_attribute_match(element, attributes):
    """Check an element to see if it matches the given attributes.

    If an element passes the check, give it a "quality" corresponding to
    the number of matched attributes. Otherwise, return null match.
    """
    match_quality = 0
    for key in element.keys():
        if key in attributes and attributes[key] == element.get(key):
            match_quality += 1
        else:
            return ElementMatch()
    return ElementMatch(element, match_quality)

def best_match(xml_tree, path, attributes={}):
    """Find the best match for a path with attributes in an XML tree.

    The return value is the match's text.

    Arguments:
    xml_tree - A tree from the xml.etree.ElementTree module.
    path - The path to search for.
    attributes - A dictionary containing attributes to match. Not all
                 attributes must be present on an element to get a match,
                 but any attributes present must match this input, and
                 the "best" match has the most matches.

    The attributes argument defaults to an empty dict.
    """

    # Recursive function to find a match.
    def find_best_below(element, path):
        # Done when there's no more of the path to match.
        if len(path) == 0:
            return ElementMatch(element)

        # Otherwise, get the beginning part of the path.
        path_head, slash, path_tail = path.partition("/")

        # Search through subelements that match the next part of the path.
        # For each one, get the quality of the match based on the number
        # of matching attributes, then call self recursively to get
        # best match from subelements.
        best_match = ElementMatch()
        for trial_element in element.findall(path_head):
            local_match = element_attribute_match(trial_element, attributes)
            if local_match:
                new_match = \
                    find_best_below(trial_element, path_tail)
                new_match += local_match
                if best_match < new_match:
                    best_match = new_match

        return best_match

    # Do the search.
    element = xml_tree.getroot()
    match = find_best_below(element, path)

    return match.get_text()

def all_matches(xml_tree, path, attributes={}):
    """Find all matches for a path with attributes in an XML tree.

    This is a generator. Each of the returned values will be the text of
    one of the matches found.

    Arguments:
    xml_tree - A tree from the xml.etree.ElementTree module.
    path - The path to search for.
    attributes - A dictionary containing attributes to match. Not all
                 attributes must be present on an element to get a match,
                 but any attributes present must match this input.

    The attributes argument defaults to an empty dict.
    """

    # Recursive function to find all matches.
    def find_matches_below(element, path):
        # Done when there's no more of the path to match.
        if len(path) == 0:
            yield ElementMatch(element)
            return

        # Otherwise, get the beginning part of the path.
        path_head, slash, path_tail = path.partition("/")

        # Search through subelements that match the next part of the path.
        # For each one, call self recursively to get all subelement matches.
        for trial_element in element.findall(path_head):
            if element_attribute_match(trial_element, attributes):
                for match in find_matches_below(trial_element, path_tail):
                    yield match

    element = xml_tree.getroot()
    for match in find_matches_below(element, path):
        yield match.get_text()
