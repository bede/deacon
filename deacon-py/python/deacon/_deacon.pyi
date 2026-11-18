from _typeshed import Incomplete
from os import PathLike
from typing import SupportsFloat, SupportsIndex, final

@final
class Index:
    """
    A loaded minimizer index, reusable across many `filter` calls.
    """
    def __new__(cls, path: str | PathLike[str], /, *, complexity_threshold: SupportsFloat | SupportsIndex | None = None) -> Index: ...
    @staticmethod
    def fetch(*, name: str = "panhuman-1", k: SupportsIndex = 31, w: SupportsIndex = 15, output: str | PathLike[str] | None = None, force: bool = False, complexity_threshold: SupportsFloat | SupportsIndex | None = None) -> Index:
        """
        Download a prebuilt index unless already present (or `force`), then load and return it.
        """
    def filter(self, input: str | PathLike[str], /, *, input2: str | PathLike[str] | None = None, interleaved: bool = False, check_pairs: bool = False, deplete: bool = False, rename: bool = False, output: str | PathLike[str] | None = None, output2: str | PathLike[str] | None = None, summary: str | PathLike[str] | None = None, abs_threshold: SupportsIndex = 2, rel_threshold: SupportsFloat | SupportsIndex = 0.01, prefix_length: SupportsIndex = 0, discard_quality: bool = False, ordered: bool = False, threads: SupportsIndex = 8, compression_level: SupportsIndex = 2, compression_threads: SupportsIndex = 0, cbq_block_size: SupportsIndex = 16, debug: str | PathLike[str] | None = None) -> dict: ...
    def info(self, /) -> dict:
        """
        Index metadata: k, w, format and minimizer/key count.
        """

def __getattr__(name: str) -> Incomplete: ...
