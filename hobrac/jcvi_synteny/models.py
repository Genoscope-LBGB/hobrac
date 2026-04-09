from dataclasses import dataclass

DEFAULT_COLOR = "lightgrey"


@dataclass
class BuscoGene:
    busco_id: str
    chromosome: str
    start: int
    end: int


@dataclass
class PairwiseAssociation:
    """Raw pairwise significant association before ALG grouping."""

    species1: str
    species2: str
    chr1: str
    chr2: str
    p_value: float
    gene_count: int


@dataclass
class ALGAssociation:
    chr1: str
    chr2: str
    p_value: float
    color: str
    gene_count: int
    alg_id: int = -1


# 37-color palette for ALG visualization
ALG_PALETTE = [
    "#e6194b",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#46f0f0",
    "#f032e6",
    "#bcf60c",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#9a6324",
    "#fffac8",
    "#800000",
    "#aaffc3",
    "#808000",
    "#ffd8b1",
    "#000075",
    "#808080",
    "#ff6f61",
    "#6b5b95",
    "#88b04b",
    "#f7cac9",
    "#92a8d1",
    "#955251",
    "#b565a7",
    "#009b77",
    "#dd4124",
    "#d65076",
    "#45b8ac",
    "#efc050",
    "#5b5ea6",
    "#9b2335",
    "#dfcfbe",
    "#55b4b0",
    "#e15d44",
]
