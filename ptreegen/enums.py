def enum(**enums):
    return type('Enum', (), enums)

TreeBuildAlgorithms = enum(
    NJ="Neigbor-joining"
    , PARSIMONY="Parsimony"
)

SeqTypes = enum(
    DNA="DNA sequences"
    , AA="Amino acid sequences"
)


