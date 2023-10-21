# -*- coding: utf-8 -*-

"""Methods to support optional packages."""

import importlib

try:
    import pytest
except ImportError:
    have_pytest = False
else:
    have_pytest = True

from whitematteranalysis.utils.tripwire import TripWire


def optional_package(name, trip_msg=None):
    """Return package-like thing and module setup for package ``name``.

    Parameters
    ----------
    name : str
        Package name.
    trip_msg : None or str
        Message to be shown when the specified package cannot be imported.

    Returns
    -------
    pckg : package, module or ``TripWire`` instance
        If the package can be imported, return it. Otherwise, return an object
        raising an error when accessed.
    have_pkg : bool
        True if import for package was successful, False otherwise.
    module_setup : function
        Callable usually set as ``setup_module`` in calling namespace, to allow
        skipping tests.

    Examples
    --------
    Typical use:

    >>> from whitematteranalysis.utils.opt_pckg import optional_package
    >>> pckg, have_pckg, _setup_module = optional_package("not_a_package")

   In this case the package doesn't exist, and therefore:

    >>> have_pckg
    False

    and

    >>> pckg.some_function() #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    TripWireError: We need package not_a_package for these functions, but
    ``import not_a_package`` raised an ImportError

    If the module does exist - we get the module

    >>> pckg, _, _ = optional_package("os")
    >>> hasattr(pckg, "path")
    True

    Or a submodule if that's what we asked for

    >>> subpckg, _, _ = optional_package("os.path")
    >>> hasattr(subpckg, "dirname")
    True
    """

    try:
        pckg = importlib.import_module(name)
    except ImportError:
        pass
    else:  # import worked
        # top level module
        return pckg, True, lambda: None
    if trip_msg is None:
        trip_msg = (
            f"{name} needed, but ``import f{name}`` raised an ``ImportError``.")
    pckg = TripWire(trip_msg)

    def setup_module():
        if have_pytest:
            pytest.mark.skip(f"No {name} for these tests.")

    return pckg, False, setup_module
