from _typeshed import Incomplete
from os import PathLike
from typing import final

@final
class Index:
    """
    A loaded minimizer index, reusable across many `filter` calls.
    """
    def __new__(cls, path: str |PathLike[str], /, *, complexity_threshold: float |None = None) -> Index: ...
    @staticmethod
    def fetch(*, name: str = "panhuman-1", k: int = 31, w: int = 15, output: str |PathLike[str] |None = None, complexity_threshold: float |None = None) -> Index:
        """
        Download a prebuilt index, then load and return it.
        """
    def filter(self, input: str |PathLike[str], /, *, input2: str |PathLike[str] |None = None, interleaved: bool = False, check_pairs: bool = False, deplete: bool = False, rename: bool = False, output: str |PathLike[str] |None = None, output2: str |PathLike[str] |None = None, summary: str |PathLike[str] |None = None, abs_threshold: int = 2, rel_threshold: float = 0.01, prefix_length: int = 0, discard_quality: bool = False, ordered: bool = False, threads: int = 8, compression_level: int = 2, compression_threads: int = 0, cbq_block_size: int = 16, quiet: bool = True, debug: bool = False) -> dict: ...
    def info(self, /) -> dict:
        """
        Index metadata: k, w, format and minimizer/key count.
        """

def __getattr__(name: str) -> Incomplete: ...
