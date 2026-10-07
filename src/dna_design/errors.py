"""Explicit domain errors; never silently discard biological information."""


class DNAError(ValueError):
    """Invalid input or unsupported biological interpretation."""


class ParseError(DNAError):
    """A construct cannot be imported faithfully."""


class UnsupportedLocationError(ParseError):
    """The simple-interval model cannot represent this location."""


class TranslationError(DNAError):
    """A CDS cannot be translated under its explicit interpretation."""


class InvalidDesignError(DNAError):
    """Final validation failed: this design must not be exported as valid."""
