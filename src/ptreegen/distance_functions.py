## @package distance_functions
# Contains implementations of possible distance measures.
#

import math

from enums import SeqTypes

##
# Computation of the uncorrected p-distance.
#
# The formula is similar to the one used in EMBOSS
# (see http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/distmat.html).
#
# @param seq1 first sequence
# @param seq2 second sequence
# @param *args positional arguments
# @param **kwargs keyword arguments
# (the "gap_penalty" argument is used to determine the gap penalty)
#
# @return distance as a single value
def p_distance(seq1, seq2, *args, **kwargs):
    assert len(seq1) == len(seq2)
    gap_penalty = 0
    if kwargs.has_key("gap_penalty"):
        gap_penalty = kwargs["gap_penalty"]
    positions_all = len(seq1)
    matches = 0
    gaps = 0
    for idx in range(positions_all):
        res1 = seq1[idx]
        res2 = seq2[idx]
        if res1 == res2 and not(res1 == "-" and res2 == "-"):
            matches+=1
        elif res1 == "-" or res2 == "-" and not(res1 == "-" and res2 == "-"):
            gaps+=1
    return 1 - float(matches) / ((positions_all - gaps) + gaps * gap_penalty)

##
# Poisson corrected p-distance.
#
# For more info see: http://goo.gl/upr3wR
#
# @param seq1 first sequence
# @param seq2 second sequence
# @param *args positional arguments
# @param **kwargs keyword arguments
#
# @return distance as a single value
def poisson_corrected(seq1, seq2, *args, **kwargs):
    p_dist = p_distance(seq1, seq2, *args, **kwargs)
    return -1 * math.log(1 - p_dist)

##
# Distance according to the Jukes-Cantor model.
#
# For more info see: http://goo.gl/upr3wR
#
# @param seq1 first sequence
# @param seq2 second sequence
# @param *args positional arguments
# @param **kwargs keyword arguments
# ("sequence_type" option is used to determine the parameters for the formula)
#
# @return distance as a single value
def jukes_cantor(seq1, seq2, *args, **kwargs):
    p_dist = p_distance(seq1, seq2, *args, **kwargs)
    if kwargs["sequence_type"] == SeqTypes.AA:
        return (-19.0/20.0) * math.log(1 - (20.0/19.0) * p_dist)
    elif kwargs["sequence_type"] == SeqTypes.AA:
        return (-3.0/4.0) * math.log(1 - (4.0/3.0) * p_dist)
    else:
        assert False
