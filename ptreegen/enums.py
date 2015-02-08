def enum(**enums):
    return type('Enum', (), enums)

TreeBuildAlgorithms = enum(
    NJ="NJ"
    , PARSIMONY="PARSIMONY"
)

SeqTypes = enum(
    DNA="DNA"
    , RNA="RNA"
    , AA="AA"
)


