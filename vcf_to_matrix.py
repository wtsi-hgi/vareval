import vcf
from typing import *
from collections import namedtuple
import itertools
from contextlib import ExitStack
import argparse
from unittest.mock import patch, Mock
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def varient_to_str(varient):
    return f"{varient.CHROM}:{varient.start + 1}-{varient.end + 1}"

def is_missing_varient_in_vcf(varient: vcf.model._Record, source_vcf_reader: vcf.Reader):
    """
    Sees if a given varient (reported as missing is the evaluation vcf) is present in a given source vcf.
    This also throws an error if the given varient is called in the source vcf.
    """
    looked_up_varients = list(source_vcf_reader.fetch(varient.CHROM, varient.start, varient.end))

    if len(looked_up_varients) == 0:
        return False
    elif len(looked_up_varients) == 1:
        if looked_up_varients[0].num_unknown != 0:
            raise Exception(f"Varient {varient_to_str(varient)}"
                f" is not called in the evaluation vcf, but is in {source_vcf_reader.filename}")
        else:
            return True
    else:
        raise Exception(f"More than one varient found with position {varient_to_str(varient)} in {source_vcf_reader.filename}")

def get_dosage_matrix_entry(varients, truth_vcf, query_vcf):
    pass

def get_dosage_matrix(eval_vcf: str, truth_vcf: str, query_vcf: str) -> List[List[int]]:
    """
    Creates 4*4 matrix with the query vcf dosages as rows and truth vcf as columns (0 indexed).
    The last column and row also represents non called varients.

    NOTE: the eval vcf must be the result of `rtg vcfeval` using the output mode of ga4gh.
    """

    matrix = [[0]*4 for _ in range(4)]

    vcf_eval_dosage = 0
    set_truth_zero = 0
    many_reprs = 0

    eval_vcf_reader = vcf.Reader(filename=eval_vcf)
    truth_vcf_reader = vcf.Reader(filename=truth_vcf)
    query_vcf_reader = vcf.Reader(filename=query_vcf)

    for variants in chunk_by(eval_vcf_reader, lambda x: x.POS):
        def increment_matrix(truth_dosage, query_dosage):
            """
            Increments the result matrix.
            NOTE: `None` represents no varient recorded
            """
            logger.debug(f"Setting truth={truth_dosage} query={query_dosage} for varients {list(map(varient_to_str, variants))}")
            if truth_dosage is None:
                truth_dosage = 3
            if query_dosage is None:
                query_dosage = 3

            # Matrix has rows of test dosages and columns of truth dosages
            matrix[query_dosage][truth_dosage] += 1

        assert len(variants) != 0

        if len(variants) == 1:
            variant = variants[0]
            truth = variant.genotype("TRUTH")
            truth_eval = truth.data.BD
            truth_dosage = truth.gt_type

            query = variant.genotype("QUERY")
            query_eval = query.data.BD
            query_dosage = query.gt_type

            if truth_eval is not None and query_eval is not None:
                vcf_eval_dosage += 1
                increment_matrix(truth_dosage, query_dosage)
            elif truth_eval is None and query_eval is not None:
                if is_missing_varient_in_vcf(variant, truth_vcf_reader):
                    increment_matrix(None, query_dosage)
                else:
                    logger.debug(f"Ignoring varient {varient_to_str(variant)} (not called in truth vcf)")
            elif truth_eval is not None and query_eval is None:
                if is_missing_varient_in_vcf(variant, query_vcf_reader):
                    increment_matrix(truth_dosage, None)
                else:
                    set_truth_zero += 1
                    increment_matrix(truth_dosage, 0)
            else:
                logger.warning(f"Ignoring varient {varient_to_str(variant)} (truth and query missing in eval vcf)")
        elif len(variants) == 2:
            bds = list(map(
                lambda x: variants[x[0]].genotype(x[1]).data.BD,
                itertools.product((1, 0), ("TRUTH", "QUERY"))
            ))
            gt_types = list(map(
                lambda x: variants[x[0]].genotype(x[1]).gt_type,
                itertools.product((1, 0), ("TRUTH", "QUERY"))
            ))

            if not (
                (bds[0], bds[3]) == (None, None) and None not in (bds[1], bds[2]) or
                (bds[1], bds[2]) == (None, None) and None not in (bds[0], bds[3])
                ):
                raise Exception(f"Invalid double varient record at {varient_to_str(variants[0])} and {varient_to_str(variants[1])}")

            truth_dosage = gt_types[0] if gt_types[2] is None else gt_types[2]
            query_dosage = gt_types[1] if gt_types[3] is None else gt_types[3]

            many_reprs += 1
            increment_matrix(truth_dosage, query_dosage)
        else:
            raise Exception(f"Many varients with the same postions. First one is at {varient_to_str(variants[0])}")

    print(f"vcf_eval_dosage: {vcf_eval_dosage}, set_truth_zero: {set_truth_zero}, many_reprs: {many_reprs}")

    return matrix


x = TypeVar("x")
def chunk_by(iterable: Iterable[x], predicate: Callable[[x], Any]) -> Iterable[List[x]]:
    iterator = iter(iterable)
    chunk = [next(iterator)]
    chunk_property = predicate(chunk[0])

    for entry in iterator:
        current_property = predicate(entry)

        if chunk_property == current_property:
            chunk.append(entry)
        else:
            yield chunk
            chunk = [entry]
            chunk_property = current_property

    yield chunk

def test_chunk_by():
    is_even = lambda x: x % 2 == 0
    assert list(chunk_by(iter([2, 4, 5, 2, 8]), is_even)) == [[2, 4], [5], [2, 8]]

    assert list(chunk_by(iter([]), is_even)) == []

def main(eval_vcf_name: str, truth_vcf_name: str, query_vcf_name: str):
    matrix = get_dosage_matrix(eval_vcf_name, truth_vcf_name, query_vcf_name)

    line_lables = ["0", "1", "2", "."]
    print(" ", *line_lables, sep="\t")
    print(*(["-"] * 5), sep="\t")
    for i, row in enumerate(matrix):
        print(i, *row, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--eval_vcf", "-e", required=True)
    parser.add_argument("--truth_vcf", "-t", required=True)
    parser.add_argument("--query_vcf", "-q", required=True)
    parser.add_argument("--debug", action="store_true")
    opts = parser.parse_args()

    if opts.debug:
        logger.setLevel(logging.DEBUG)

    main(opts.eval_vcf, opts.truth_vcf, opts.query_vcf)
