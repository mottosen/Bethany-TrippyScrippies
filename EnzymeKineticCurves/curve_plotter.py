import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tools.michaelis_menten import MichaelisMenten
from tools.utility import Unit, UtilityUnit

cur_dir = Path(".")


class Config:
    """This class handles parsing of cmd input."""

    def __init__(self):
        self._args = sys.argv
        self._data_dir = cur_dir / "data"
        self._res_dir = cur_dir / "results"
        self._unit = UtilityUnit.Micro
        self._extinction_coefficient = 0.021
        self._time_coefficient = 95 / 60
        self._save_fig = False
        self._plotter = None
        self._title = "Title Not Specified"

        self._parseArguments(sys.argv)

        if not self._res_dir.exists():
            os.makedirs(self._res_dir)

    def GetDataDir(self):
        return self._data_dir

    def GetResultsDir(self):
        return self._res_dir

    def GetUnit(self):
        return self._unit

    def GetExtinctionCoefficient(self):
        return self._extinction_coefficient

    def GetTimeCoefficient(self):
        return self._time_coefficient

    def GetSaveMode(self):
        return self._save_fig

    def GetPlotter(self):
        return self._plotter

    def GetPlotTitle(self):
        return self._title

    @staticmethod
    def GetHelpText():
        return """
Usage: ([..] is optional)
    curve_plotter.py [options] action

    action:
        MichaelisMenten

    options:
        -unit <unit specifier>  => the unit used, default=micro
        -title <plot title>     => the title of plot being created
        -save                   => whether to save the plot, default=show

    unit specifier:
        milli, micro



    example use:
        curve_plotter.py -save -unit "milli" -title "EnzymeXYZ plotted for DATA"

        this will save a created plot, with the specified title, using milli moles as unit.
"""

    def _parseArguments(self, args):
        mode = args[-1]
        options = args[1:-1]
        i = 0
        lim = len(options)

        while i < lim:
            match options[i]:
                case "-save":
                    self._save_fig = True
                    i += 1

                case "-unit":
                    self._unit = UtilityUnit.from_text(options[i + 1])
                    i += 2

                case "-title":
                    self._title = options[i + 1]
                    i += 2

                case _:
                    raise Exception("Could not decide on option.")

        match mode:
            case "MichaelisMenten":
                self._plotter = MichaelisMenten(self._unit)

            case _:
                raise Exception("Could not decide on plotter.")


class DataReader:
    def __init__(self, config):
        self._config: Config = config

    def _handle_column(self, col: [float]) -> Unit:
        """
        handle_column

        Handles a column of a data file.

        First, normalizes the column. Then, computes a linear fit.
        """

        # normalize the column
        init_val = col[0]
        col_norm = [entry - init_val for entry in col]

        # linear fit for data points
        slope, _ = np.polyfit(list(range(len(col_norm))), col_norm, 1)

        # reduction of substrate per minute
        init_velocity = (
            slope / config.GetTimeCoefficient() / config.GetExtinctionCoefficient()
        )

        return init_velocity

    def handle_file(self, file_path: str, config: Config) -> (Unit, [float]):
        """
        handle_file

        Handles a data file.
        """

        # extract concentration from file name
        concentration = float(
            re.search(r"(\d+(_(\d+))?)", str(file_path)).group(0).replace("_", ".")
        )

        # read data and compute initial velocity
        df = pd.read_csv(file_path, header=None)
        reduction_cytochromeC = df.apply(self._handle_column, axis=0)

        return Unit(concentration), reduction_cytochromeC


if __name__ == "__main__":
    if sys.argv[1] in ["-h", "--help"]:
        print(Config.GetHelpText())
        exit()

    config = Config()
    reader = DataReader(config)

    data = {}

    for data_file in config.GetDataDir().glob("*"):
        concentration, velocity = reader.handle_file(data_file, config)

        data[concentration.get_num(config.GetUnit())] = velocity

    data = dict(sorted(data.items()))

    plotter = config.GetPlotter()
    plotter.make_plot(
        data, config.GetResultsDir(), config.GetPlotTitle(), config.GetSaveMode()
    )
