## @package enums
# Contains all enums used in the project.
#

##
# Definition of the enum "type".
#
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

OutputForm = enum(
    NEWICK = "NEWICK",
    IMAGE_PNG = "IMAGE_PNG",
    IMAGE_SVG = "IMAGE_SVG",
    GUI = "GUI",
    PRINT = "PRINT",
)

TreeType = enum(
    CIRC = "CIRC",
    RECT = "RECT",
)

##
# @var TreeBuildAlgorithms
# This enum contains possible algorithm choices.
# @var SeqTypes
# This enum contains possible sequence types.
# @var DistMeasures
# This enum contains possible distance measures.
#


