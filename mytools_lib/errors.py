"""Exceptions shared by mytools subcommands.

Caught at the top level (mytools) and turned into the exit codes from
SPEC.md sec7: BedError -> 1 (bad input data), UsageError -> 2 (usage error).
"""


class BedError(Exception):
    pass


class UsageError(Exception):
    pass
