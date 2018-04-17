from typing import Iterable, IO, Any
from collections import Counter, namedtuple
from itertools import chain

test_data = \
"""10 01 N0
00 10 01
00 00 00"""

test_data2 = \
"""0.123 10 01 N0
0.12 00 10 01
0.14 00 00 00"""

Stats = namedtuple("Stats", [
    "concordance",
    "discordance",
    "x_na",
    "y_na"
])

def get_stats(data: Iterable[str]):
    counter = Counter(data)
    concordance = 0
    discordance = 0
    x_na = 0
    y_na = 0

    for key in counter.keys():
        assert len(key) == 2
        if key[0] == key[1]:
            concordance += counter[key]
        elif "N" in key[0]:
            x_na += counter[key]
        elif "N" in key[1]:
            y_na += counter[key]
        else:
            discordance += counter[key]

    return Stats(
        concordance,
        discordance,
        x_na,
        y_na
    )

def test_get_stats():
    stats = get_stats(["01", "01", "11", "N1"])

    assert stats.concordance == 1
    assert stats.discordance == 2
    assert stats.x_na == 1
    assert stats.y_na == 0

def flattern(lst: Iterable[Iterable[Any]]):
    return chain(*lst)

def display_stats(text: str, interval_size: int):
    num_intervals = int(100/interval_size)
    assert 100/num_intervals == interval_size

    intervals_results = [[] for x in range(num_intervals)]
    for line in text.split("\n"):
        columns = line.split(" ")
        allele_freq = float(columns[0])
        assert 0 <= allele_freq <= 1
        if allele_freq == 1:
            intervals_results[-1].extend(columns[1:])
        else:
            print(intervals_results)
            intervals_results[int(allele_freq / interval_size)].extend(columns[1:])

    for i, result in enumerate(intervals_results):
        print(f"[{i*interval_size}-{(i+1)*interval_size}]: {get_stats(result)}")

display_stats(test_data2, 5)