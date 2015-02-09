
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

DistMeasures = enum(
    P_DISTANCE = "P_DISTANCE",
    POISSON_CORRECTED = "POISSON_CORRECTED",
    JUKES_CANTOR = "JUKES_CANTOR",
)


