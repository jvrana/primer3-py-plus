class Primer3PlusException(Exception):
    """Generic primer3plus exception."""


class Primer3PlusParserError(Primer3PlusException):
    """Generic parser exception."""


class Primer3PlusRunTimeError(Exception):
    """Exception for errors returned from primer3."""


class Primer3PlusWarning(Warning):
    """Warning for primer3plus."""
