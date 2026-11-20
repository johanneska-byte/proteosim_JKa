import proteosim_JKa as ps
from pyteomics import achrom
def test_predict_lc_retention_times():
    peptides = [
        "HCILNKLTEWL",
        "AAAAAARRRRRGHHK",
        "RRRRRAW"
    ]

    expected = {
        "HCILNKLTEWL": 47.6,
        "AAAAAARRRRRGHHK": 21.70,
        "RRRRRAW": 16.20
    }

    actual = ps.predict_lc_retention_times(peptides)

    actual = {
        peptide: round(rt, 2)
        for peptide, rt in actual.items()
    }

    assert actual == expected

def test_select_retention_time_window():

    peptide_rt_map = {
        "PEPTIDE": 12.3,
        "AAAKR": 18.7,
        "RRRRRAW": 25.1
    }

    selected = ps.select_retention_time_window(
        peptide_rt_map,
        lower_ret_time=15,
        upper_ret_time=25
    )

    window_dict_dummy = {}

    for seq, rt in peptide_rt_map.items():
        if 15 <= rt <= 25:
            window_dict_dummy[seq] = rt

    assert selected == window_dict_dummy
