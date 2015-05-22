
import pandas as pd


try:
    from functools import lru_cache
except ImportError:
    # dummy cache, i.e. Python 2 version will be a bit slower.
    def lru_cache():
        def dummy(func):
            return func

        return dummy


class AbstractResults:
    def __init__(self, dataframe):
        """
        Arguments:
        dataframe -- a pandas data frame or its path consisting of per gene results as produced by MAGeCK
        """
        self.df = pd.read_table(dataframe, na_filter=False)

    def __getitem__(self, slice):
        return self.df.__getitem__(slice)
