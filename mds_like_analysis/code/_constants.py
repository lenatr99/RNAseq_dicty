"""
This file contains all the constants used in the analysis.
"""

DATA_PATH_ALL = "../data/expression_data.csv"
DATA_PATH = "../data/expression_data_average.csv"
ANNOT_PATH = "../data/DictyGeneAnnotations_3504.csv"
SCALING = ["None", "m0s1"]
NAME_DICT = {
    "AX4": "$AX4$",
    "B1-": "$tgrB1^-$",
    "C1-": "$tgrC1^-$",
    "rgB-": "$rapgapB^-$",
    "B1-rgB-": "$rapgapB^-tgrB1^-$",
    "AX4L846F": "$AX4$ $L846F$",
}
LINE_COLORS = {
    "AX4": "#000000",
    "B1-": "#00B2FF",
    "C1-": "#A400D3",
    "rgB-": "#008528",
    "B1-rgB-": "#D9D800",
    "AX4L846F": "#ED1C24",
}

MILESTONES = [
    "noagg_ripple",
    "ripple_lag",
    "lag_tag",
    "tag_tip",
    "tip_slug",
    "slug_Mhat",
    "Mhat_cul",
    "cul_FB",
]

MILESTONE_COLORS = {
    "noagg_ripple_up": "#FF69B4",
    "noagg_ripple_down": "#FF69B4",
    "ripple_lag_up": "#00B2FF",
    "ripple_lag_down": "#00B2FF",
    "lag_tag_up": "#A400D3",
    "lag_tag_down": "#A400D3",
    "tag_tip_up": "#FF69B4",
    "tag_tip_down": "#008528",
    "tip_slug_up": "#D9D800",
    "tip_slug_down": "#D9D800",
    "slug_Mhat_up": "#ED1C24",
    "slug_Mhat_down": "#4B0082",
    "Mhat_cul_up": "#FF7F00",
    "Mhat_cul_down": "#FF7F00",
    "cul_FB_up": "#4B0082",
    "cul_FB_down": "#4B0082",
}


SPACE = 100
SPACE_VAR = 0.5
LENGTH_DICT = {
    "AX4": 1.5,
    "B1-": 2.5,
    "C1-": 2.5,
    "rgB-": 3,
    "B1-rgB-": 5.5,
    "AX4L846F": 3.5,
}
