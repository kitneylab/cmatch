import datetime
import logging
import time
from functools import wraps
from os import path
from string import whitespace
from typing import Any, Callable
import json

from Bio import Seq, SeqRecord


def read_json(json_file):
    """
    Read JSON file

    Args:
        json_file (str): The path to the JSON file

    Returns:
        dict: JSON object

    """
    with open(json_file) as json_file:
        json_dict = json.load(json_file)
    return json_dict


def timeit(func: Callable[..., Any]) -> Callable[..., Any]:
    """Times a function, usually used as decorator"""

    @wraps(func)
    def timed_func(*args: Any, **kwargs: Any) -> Any:
        """Returns the timed function"""
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = datetime.timedelta(seconds=(time.time() - start_time))
        logging.info(f"Time spent on {func.__name__}(): {elapsed_time}")
        return result

    return timed_func


def biopython_stfu():
    """
    Ignore biopython warnings
    """
    import warnings
    from Bio import BiopythonWarning

    warnings.simplefilter("ignore", BiopythonWarning)


def read_file(filename: str):
    """
    Read file and return content

    :param filename:

    """
    with open(filename, "r") as file_object:
        file_contents = file_object.read()
    return file_contents


def file_to_seqrec(filename: str):
    """
    Read file return SeqRec

    :param filename:

    """
    basename = path.basename(filename)
    name = path.splitext(basename)[0]
    seq = Seq.Seq(read_file(filename).strip(whitespace))  # strip newlines and shit
    seqrec = SeqRecord.SeqRecord(seq, id=name, name=name, description=name)
    return seqrec
