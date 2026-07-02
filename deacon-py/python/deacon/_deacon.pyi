from _typeshed import Incomplete
from typing import Any, final

@final
class Index:
    """
    A loaded minimizer index, reusable across many `filter` calls.
    """
    def __new__(cls, /, path: str, complexity_threshold: float |None = None) -> Index: ...
    @staticmethod
    def fetch(name: str = "panhuman-1", k: int = 31, w: int = 15, output: str |None = None, complexity_threshold: float |None = None) -> Index:
        """
        Download a prebuilt index, then load and return it.
        """
    def filter(self, /, fastq: str, fastq2: str |None = None, interleaved: bool = False, deplete: bool = False, rename: bool = False, rename_random: bool = False, output: str |None = None, output2: str |None = None, abs_threshold: int = 2, rel_threshold: float = 0.01, prefix_length: int = 0, output_fasta: bool = False, threads: int = 8, compression_level: int = 2, compression_threads: int = 0, debug: bool = False, quiet: bool = True) -> Any: ...
    def info(self, /) -> dict:
        """
        Index metadata: k, w, format and minimizer/key count.
        """

def __getattr__(name: str) -> Incomplete: ...
