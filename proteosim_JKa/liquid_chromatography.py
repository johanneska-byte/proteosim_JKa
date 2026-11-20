
def predict_lc_retention_times(peptides):
    """
    Takes a list of peptides and returnes a dictionary mapping each peptide to its predicted retention time.

    Parameters
    ----------
    peptides: list
        List of peptide sequences
    Returns
    -------
    dict
    """
    retention_time = {}
    for peptide in peptides:
        value = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        retention_time[peptide] = round(value, 3)

    return retention_time

from pyteomics import achrom
import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Creates a histogram of the retention times of peptides in the retention_times variable with'list' number of histogram bins. Accepts list and dict. 

    Parameters
    ----------
    retention_times : list or dict
        Retention times to plot. If a dictionary is provided, only the dictionary values are used

    resolution : int, default=30
        Number of histogram bins

    Returns
    -------
    None
        Displays a histogram of the retention time distribution
    
    """
    plt.hist(retention_times, bins = resolution, color= "#666666")
    plt.xlabel("Retention time")
    plt.ylabel("count")
    plt.title("Retention time distribution")

    plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Filters a peptide retention time dictionary for peptides whose predicted retention times fall within a specified retention time window

    Parameters
    ----------
    peptide_rt_map : dict
        Dictionary mapping peptide sequences to their predicted retention times

    lower_ret_time : float
        Lower bound of the retention time window.

    upper_ret_time : float
        Upper bound of the retention time window

    Returns
    -------
    dict
        Dictionary containing only peptide/retention-time pairs with retention times between lower_ret_time and upper_ret_time
    """
    window_dict_gatekept = {}
    for seq, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time: #gatekeeping
            window_dict_gatekept[seq] = rt

    return window_dict_gatekept