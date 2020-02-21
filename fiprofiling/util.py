#!/usr/bin/env python3
import sys
import itertools

import pandas as pd

def bool2binstring(booliter):
    if isinstance(booliter, bool):
        return '1' if booliter else '0'
    binstring = ''.join(('1' if i else '0' for i in booliter))
    return binstring

def binstring2bool(boolstring):
    return tuple(i == '1' for i in boolstring)

def stdin2binstrings():
    binstring = next(sys.stdin)
    # check that the first one consists only of 0's and 1's
    binstring = binstring.strip()
    if not set(binstring) <= set(('0', '1')):
        raise ValueError(
            "The provided string '{0}' contains different symbols than just 0 and 1.".format(
                binstring))
    yield binstring
    for further_binstring in sys.stdin:
        yield further_binstring.strip()

def stdin2dataframe():
    return pd.DataFrame.from_csv(sys.stdin)

def chunk_iterable(iterable, chunksize=1000):
    iterator = iter(iterable)
    chunk = list(itertools.islice(iterator, chunksize))
    while chunk:
        yield chunk
        chunk = list(itertools.islice(iterator, chunksize))
