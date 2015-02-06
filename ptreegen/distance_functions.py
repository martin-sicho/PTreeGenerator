
def p_distance(seq1, seq2, gapPenalty=0):
    assert len(seq1) == len(seq2)
    positions_all = len(seq1)
    matches = 0
    gaps = 0
    for idx in range(positions_all):
        res1 = seq1[idx]
        res2 = seq2[idx]
        if res1 == res2 and not(res1 == "-" and res2 == "-"):
            matches+=1
        elif res1 == "-" or res2 == "-":
            gaps+=1
    return 1 - float(matches) / ((positions_all - gaps) + gaps * gapPenalty)