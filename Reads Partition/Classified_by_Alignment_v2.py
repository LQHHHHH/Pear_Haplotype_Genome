#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classify reads between two BAMs based on alignment score (AS tag).

For each read present in BOTH BAMs:
- If mean(AS in bam1) > mean(AS in bam2): assign to hap1 only
- If mean(AS in bam1) < mean(AS in bam2): assign to hap2 only
- If equal: assign to BOTH hap1 and hap2

Works with BAMs in any sort order.

Outputs:
  <prefix1>.hap1.id  - reads assigned to hap1 (including ties)
  <prefix2>.hap2.id  - reads assigned to hap2 (including ties)

Printed summary includes:
  - total hap1/hap2 reads (including ties)
  - numbers of hap1-only / hap2-only reads
  - ratios between them
"""

import argparse
import logging
from statistics import mean
from typing import Dict, List, Tuple

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s"
)


class ReadBam:
    """Read a BAM and collect AS scores per read name."""

    def __init__(self, bam: pysam.AlignmentFile) -> None:
        self.bam = bam

    def get_paired_scores(self) -> Dict[str, List[int]]:
        """
        Return mapping: read_name -> [AS1, AS2, ...]

        Only records that have an 'AS' tag are used.
        """
        scores: Dict[str, List[int]] = {}
        missing_as = 0

        for aln in self.bam:
            try:
                as_score = aln.get_tag("AS")
            except KeyError:
                missing_as += 1
                continue
            scores.setdefault(aln.query_name, []).append(int(as_score))

        if missing_as:
            logging.warning(
                "Skipped %d alignments without AS tag in %s",
                missing_as,
                getattr(self.bam, "filename", "BAM")
            )

        return scores


class CompareBams:
    """
    Compare two BAM-derived score dictionaries.

    bam1_dict / bam2_dict are from ReadBam.get_paired_scores()
    """

    def __init__(self, bam1_dict: Dict[str, List[int]], bam2_dict: Dict[str, List[int]]) -> None:
        set1, set2 = set(bam1_dict.keys()), set(bam2_dict.keys())

        self.reads_only_1 = set1 - set2
        self.reads_only_2 = set2 - set1
        self.reads_both = set1 & set2

        self.bam1_dict = bam1_dict
        self.bam2_dict = bam2_dict

    @staticmethod
    def _compare_scores(scores1: List[int], scores2: List[int]) -> str:
        """
        Compare mean scores between two lists.

        Returns:
            'a' if mean(scores1) > mean(scores2)
            'b' if mean(scores1) < mean(scores2)
            'c' if equal
        """
        m1 = mean(scores1)
        m2 = mean(scores2)

        if m1 > m2:
            return "a"
        elif m1 < m2:
            return "b"
        else:
            return "c"

    def classify(self, prefix1: str, prefix2: str) -> Tuple[int, int, int]:
        """
        Classify reads and write ID lists.

        Parameters
        ----------
        prefix1 : str
            Output prefix for hap1 ID file.
        prefix2 : str
            Output prefix for hap2 ID file.

        Returns
        -------
        (hap1_only, hap2_only, both)
        """
        hap1_only = 0
        hap2_only = 0
        both = 0

        out1_path = f"{prefix1}.hap1.id"
        out2_path = f"{prefix2}.hap2.id"

        logging.info("Writing hap1 IDs to %s", out1_path)
        logging.info("Writing hap2 IDs to %s", out2_path)

        with open(out1_path, "w") as out1, open(out2_path, "w") as out2:
            for read_name in self.reads_both:
                scores1 = self.bam1_dict[read_name]
                scores2 = self.bam2_dict[read_name]

                label = self._compare_scores(scores1, scores2)
                if label == "a":
                    hap1_only += 1
                    out1.write(read_name + "\n")
                elif label == "b":
                    hap2_only += 1
                    out2.write(read_name + "\n")
                else:  # 'c' – equal
                    both += 1
                    out1.write(read_name + "\n")
                    out2.write(read_name + "\n")

        return hap1_only, hap2_only, both


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Classify reads between two BAMs based on alignment score (AS tag)."
    )
    parser.add_argument("bam1", help="BAM file for haplotype 1.")
    parser.add_argument("bam2", help="BAM file for haplotype 2.")
    parser.add_argument(
        "--prefix1",
        default=None,
        help="Output prefix for hap1 ID file (default: bam1 filename).",
    )
    parser.add_argument(
        "--prefix2",
        default=None,
        help="Output prefix for hap2 ID file (default: bam2 filename).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    prefix1 = args.prefix1 or args.bam1
    prefix2 = args.prefix2 or args.bam2

    logging.info("Opening BAM1: %s", args.bam1)
    bam1 = pysam.AlignmentFile(args.bam1, "rb")

    logging.info("Opening BAM2: %s", args.bam2)
    bam2 = pysam.AlignmentFile(args.bam2, "rb")

    logging.info("Reading BAM1 scores...")
    dct1 = ReadBam(bam1).get_paired_scores()

    logging.info("Reading BAM2 scores...")
    dct2 = ReadBam(bam2).get_paired_scores()

    logging.info("Comparing BAMs...")
    cmp_bams = CompareBams(dct1, dct2)
    hap1_only, hap2_only, both = cmp_bams.classify(prefix1, prefix2)

    hap1_total = hap1_only + both
    hap2_total = hap2_only + both

    # Avoid ZeroDivisionError
    ratio_total = float("inf") if hap2_total == 0 else hap1_total / hap2_total
    ratio_specific = float("inf") if hap2_only == 0 else hap1_only / hap2_only

    print(
        (
            "\n"
            f"Hap1 reads (including ties): {hap1_total}\n"
            f"Hap2 reads (including ties): {hap2_total}\n"
            f"Common (assigned to both):  {both}\n"
            f"Hap1-specific reads:        {hap1_only}\n"
            f"Hap2-specific reads:        {hap2_only}\n"
            f"hap1_vs_hap2_total_ratio:   {ratio_total:.2f}\n"
            f"specific_ratio_hap1_vs_hap2:{ratio_specific:.2f}\n"
        )
    )


if __name__ == "__main__":
    main()