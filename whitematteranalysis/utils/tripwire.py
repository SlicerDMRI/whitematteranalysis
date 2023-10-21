# -*- coding: utf-8 -*-

class TripWireError(AttributeError):
    """Class to raise an exception for missing modules or other misfortunes.."""


def is_tripwire(obj):
    """Returns True if ``obj`` appears to be a ``TripWire`` object.

    Examples
    --------
    >>> is_tripwire(object())
    False
    >>> is_tripwire(TripWire("some message"))
    True
    """

    try:
        obj.any_attribute
    except TripWireError:
        return True
    except Exception:
        pass
    return False


class TripWire:
    """Class raising error if used.

    Standard use is to proxy modules that could not be imported.

    Examples
    --------
    >>> try:
    ...     import not_a_package
    ... except ImportError:
    ...    not_a_package_name = TripWire("Do not have not_a_package_name")
    >>> not_a_package_name.do_something("with argument") #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    TripWireError: Do not have not_a_package_name
    """

    def __init__(self, msg):
        self._msg = msg

    def __getattr__(self, attr_name):
        """Raise informative error accessing attributes.
        """
        raise TripWireError(self._msg)

    def __call__(self, *args, **kwargs):
        """Raise informative error while calling.
        """
        raise TripWireError(self._msg)
