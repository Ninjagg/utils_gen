# -*- coding: utf-8 -*-
# ===============================================================================
# Data      : 20250317
# Author    : Xuanming
# Annotation: This is for creating the object "PhysicoChem" and calculating molecular_weight, isoelectric_point, and extinction_coefficient.
# ===============================================================================
import re
from typing import Union

from Bio.SeqUtils.ProtParam import ProteinAnalysis

class PhysicoChem:

    def __init__(self):
        pass

    @staticmethod
    def cal_molecular_weight(query_sequence: str) -> float:
        query_sequence = query_sequence.replace("X", "")
        molecular_weight = ProteinAnalysis(query_sequence).molecular_weight() / 1000
        return molecular_weight

    @staticmethod
    def cal_isoelectric_point(query_sequence: str) -> float:
        query_sequence = re.sub(r"[BJOUXZ]|\s", "", query_sequence)
        pI = ProteinAnalysis(query_sequence).isoelectric_point()
        return pI

    @staticmethod
    def cal_extinction_coefficient(query_sequence: str, molecular_weight: Union[int, float]) -> float:
        query_sequence = re.sub(r"[BJOUXZ]", "", query_sequence)
        amino_acid_c_number = int(query_sequence.count("C") / 2)
        amino_acid_w_number = query_sequence.count("W")
        amino_acid_y_number = query_sequence.count("Y")
        coef = (
            (
                amino_acid_w_number * 5500
                + amino_acid_y_number * 1490
                + amino_acid_c_number * 125
            )
            / molecular_weight
            / 1000
        )
        return coef