"""Explicit conversions between UI and internal coordinate conventions."""


def _bounds(start: int, end: int, length: int | None, *, ui: bool) -> None:
    if type(start) is not int or type(end) is not int:
        raise ValueError("Coordinates must be integers, not booleans or floats")
    if length is not None and (type(length) is not int or length < 1):
        raise ValueError("Sequence length must be a positive integer")
    if ui:
        valid = 1 <= start <= end
    else:
        valid = 0 <= start < end
    if not valid or (length is not None and end > length):
        raise ValueError("Coordinates are empty, reversed, or out of bounds")


def ui_to_internal(start: int, end: int, length: int | None = None) -> tuple[int, int]:
    """Convert a nonempty 1-based inclusive interval to 0-based half-open."""
    _bounds(start, end, length, ui=True)
    return start - 1, end


def internal_to_ui(start: int, end: int, length: int | None = None) -> tuple[int, int]:
    """Convert a nonempty 0-based half-open interval to 1-based inclusive."""
    _bounds(start, end, length, ui=False)
    return start + 1, end
