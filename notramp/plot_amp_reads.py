import logging
from typing import Union

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as e:
    logger.warning("Module matplotlib could not be be found.")

def gen_overview_fig(data1: dict, data2: dict, kw: dict):
    fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8, 10))
    print(data1)
    print(data2)
    plot_reads_per_amp(ax1, data1, "Available reads")
    plot_reads_per_amp(ax2, data2, "Selected reads", helpline=kw["figures"])
    fig.tight_layout()
    plt.savefig(f"{kw['out_dir']}/notramp_amplicon_coverage.png")


def plot_reads_per_amp(
        ax: object,
        data: dict,
        title: str,
        helpline: Union[bool, int]=False,
        ):
    categories = list(data.keys())
    values = list(data.values())  
    ax.bar(categories, values)
    ax.set_xlabel('Amplicon')
    ax.set_ylabel('Count')
    ax.set_title(title)
    if int(helpline):
        ax.axhline(
            y=int(helpline), color='red', linestyle='--', label=helpline
            )
