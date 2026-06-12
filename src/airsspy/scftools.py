"""
SCF convergence analysis tools for CASTEP calculations.

Provides the SCFInfo class for parsing, summarising, and plotting
self-consistent field (SCF) convergence data from .castep output files.
"""

import re
from collections import namedtuple

import numpy as np


class SCFInfo:
    """Class for storing and extracting SCF convergence information.

    Parses SCF loop data (energy, Fermi energy, energy gain, timing)
    from a CASTEP .castep file and provides summary statistics and
    plotting methods.
    """

    SCF_LINE = re.compile(
        r"^ +([0-9]+) +([+-.E0-9]+) +([+-.E0-9]+) +([+-.E0-9]+) +([.0-9]+) +<-- SCF"
    )
    ScfData = namedtuple("ScfData", ["loops", "energies", "fermi_energies", "gains", "timers"])

    def __init__(self, castep_file: str) -> None:
        """
        Construct an SCFInfo object from a .castep file.

        Args:
            castep_file: Path to the .castep file.
        """
        self.filename = castep_file
        with open(castep_file) as fhandle:
            self.scf_data: list[SCFInfo.ScfData] = self.parse(fhandle)
        self.compute_converge_data()

    def __len__(self) -> int:
        return len(self.scf_data)

    def reload(self) -> None:
        """Re-read the .castep file and reparse."""
        with open(self.filename) as fhandle:
            self.scf_data = self.parse(fhandle)
        self.compute_converge_data()

    def parse(self, lines) -> list["SCFInfo.ScfData"]:
        """
        Parse SCF loop data from lines of a .castep file.

        Args:
            lines: Iterable of lines from the .castep file.

        Returns:
            A list of ScfData namedtuples, one per geometry step.
        """
        scf_loops: list = []
        pattern = self.SCF_LINE
        loops: list[int] = []
        engs: list[float] = []
        fermi_engs: list[float] = []
        gains: list[float] = []
        timers: list[float] = []

        for line in lines:
            match = pattern.match(line)
            if not match:
                continue

            loop = int(match.group(1))
            if loop == 1:
                loops, engs, fermi_engs, gains, timers = [], [], [], [], []
                scf_loops.append(
                    self.ScfData(
                        loops=loops,
                        energies=engs,
                        fermi_energies=fermi_engs,
                        gains=gains,
                        timers=timers,
                    )
                )

            loops.append(loop)
            engs.append(float(match.group(2)))
            fermi_engs.append(float(match.group(3)))
            gains.append(float(match.group(4)))
            timers.append(float(match.group(5)))
        return scf_loops

    def compute_converge_data(self) -> None:
        """Compute summary convergence information across ionic steps."""
        energies: list[float] = []
        timings: list[float] = []
        durations: list[float] = []
        steps: list[int] = []
        avg_scf: list[float] = []
        for cycle in self.scf_data:
            init_time = cycle.timers[0]
            final_time = cycle.timers[-1]
            duration = final_time - init_time
            energies.append(cycle.energies[-1])
            timings.append(final_time)
            durations.append(duration)
            steps.append(cycle.loops[-1])
            avg_scf.append(duration / cycle.loops[-1])

        self.conv_data: dict[str, list] = {
            "energies": energies,
            "timings": timings,
            "durations": durations,
            "steps": steps,
            "average_scf_time": avg_scf,
        }

    def get_summary(self) -> dict[str, float]:
        """
        Obtain summary statistics.

        Returns:
            Dictionary with avg_ionic_time, avg_elec_time, avg_elec_steps,
            ionic_steps, and total_time.
        """
        conv_data = self.conv_data
        return {
            "avg_ionic_time": float(np.mean(conv_data["durations"])),
            "avg_elec_time": float(np.mean(conv_data["average_scf_time"])),
            "avg_elec_steps": float(np.mean(conv_data["steps"])),
            "ionic_steps": len(conv_data["steps"]),
            "total_time": sum(conv_data["durations"]),
        }

    def plot_scf(
        self, scf_no: int, xaxis: str = "loops", show: bool = True
    ):
        """
        Plot SCF convergence for a single geometry step.

        Args:
            scf_no: Index of the SCF cycle to plot.
            xaxis: Which x-axis data to use (``"loops"`` or ``"timers"``).
            show: Whether to call ``fig.show()``.

        Returns:
            The matplotlib Figure object.
        """
        import matplotlib.pyplot as plt

        data = self.scf_data[scf_no]
        fig, axs = plt.subplots(3, 1, sharex=True)
        xax = getattr(data, xaxis)
        axs[0].set_title(f"SCF Information for cycle {scf_no}")
        axs[0].plot(xax, data.energies)
        axs[0].set_ylabel("Energy (eV)")
        axs[1].plot(xax, data.fermi_energies)
        axs[1].set_ylabel("Fermi Energy (eV)")

        abs_diff = np.abs(np.array(data.gains))
        axs[2].plot(xax, abs_diff)
        axs[2].set_yscale("log")
        axs[2].set_ylabel("Energy change (eV)")
        if show:
            fig.show()
        return fig

    def plot_conv(self, axs=None, show: bool = True):
        """
        Plot convergence summary across all ionic steps.

        Args:
            axs: Optional matplotlib axes (4 subplots). Created if None.
            show: Whether to show the figure.

        Returns:
            The matplotlib Figure object (or None if axs was provided).
        """
        import matplotlib.pyplot as plt

        data = self.conv_data
        fig = None
        if not axs:
            fig, axs = plt.subplots(4, 1, sharex=True, figsize=(7, 8))
        loops = list(range(len(data["steps"])))
        axs[0].plot(loops, data["energies"], "-x", label="Final energy")
        axs[0].set_ylabel("Energy (eV)")
        axs[1].plot(loops, data["durations"], "-x", label="Durations")
        axs[1].set_ylabel("Duration (s)")
        axs[2].plot(loops, data["average_scf_time"], "-x", label="SCF time")
        axs[2].set_ylabel("Avg SCF time (s)")
        axs[3].plot(loops, data["steps"], "-x", label="Electronic steps")
        axs[3].set_ylabel("Electronic Steps")
        axs[3].set_xlabel("Ionic Step")

        if show and fig is not None:
            fig.tight_layout()
            fig.show()
        return fig
