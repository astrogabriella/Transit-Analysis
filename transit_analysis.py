#
# Gabriella Masiha
#
# For a chosen light curve from the NASA Exoplanet Archive: 
# http_seconds://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
#
# this program performs phase-folding, finds a best-fit model for the 
# phase-folded data via Chi-square testing and uses the associated parameters 
# to calculate, display and save planet parameters. You will need the 
# stellar radius[solar radii] and stellar surface gravity[log10(cm/s**2)]
# for your chosen system. 
# 
# Plot downloaded must be a Kepler DV Time Series and Periodogram 
# with LC_DETREND selected under the Y axis  column and must be in table format
#
# TBL must be in the same working directory as this script
#
# Tested on Python 3.7.9 with NumPy version 1.19.2 and Plotly version 4.12.0
# 

import os

import time
import itertools
import math
import numpy as np
import plotly.graph_objs as go
import plotly.offline as ply
from plotly.graph_objs import Figure

DEBUG_MODE = False  # Set to True while developing


def debug(*args, **kwargs):
    if DEBUG_MODE:
        print(*args, **kwargs)


def get_float(prompt: str) -> float:
    """Prompt the user for a number and converts it to a float automatically"""
    while True:
        try:
            float_ = float(input(prompt))
        except ValueError:
            print("Invalid input. Please try again.")
            continue

        return float_


def get_float_list(prompt: str, n: int) -> list:
    """Prompt the user for a list of comma-separated numbers and convert them
    to floats automatically"""

    floats = []
    while True:
        str_floats = input(prompt).split(",")
        try:
            floats = list(map(float, str_floats))
        except ValueError:
            print("Invalid input. Please try again.")
            continue
        if len(floats) != n:
            print("Invalid number of inputs. Please try again.")
            continue

        return floats


def ask_y_n(prompt) -> bool:
    """Ask a yes or no question until the user responds with a valid 
    Yes or No"""

    answer = None
    while answer is None:
        yn = input(prompt + "[y/N] ").lower()
        if yn in ("y", "yes"):
            answer = True
        elif yn in ("n", "no"):
            answer = False
        else:
            # input not recognised
            print(f"{yn!r} not valid input.")
            continue

    return answer


def make_array(start: float, end: float) -> list:
    """Return evenly spaced floats between defined limits"""
    evenly_spaced_array = np.linspace(start, end, num=11)

    return evenly_spaced_array


def make_trace(x_axis, y_axis, colour, name=""):
    """Make a trace dict for a Plotly.
    Optionally provide a name for the trace"""
    trace = {"type": "scatter",
             "mode": "markers",
             "x": x_axis,
             "y": y_axis,
             "marker": dict(
                 color=colour,
                  size=5)}
    if name:
        trace["name"] = name

    return trace

# Plotly is used for its ability to keep plots open while code is being executed
def make_plot(traces: list, kepID: int, x_axis_title: str,
              filename: str = ""):
    """Create and display plots and save them as a HTML file"""

    layout = {"title": kepID,
              "title_x": 0.5,
              "xaxis": {"title": x_axis_title},
              "yaxis": {"title": "Relative Flux"},
              "plot_bgcolor": "rgba(0,0,0,0)"
              }

    data = []
    for trace in traces:
        data.append(go.Scattergl(trace))

    fig = Figure(data=data, layout=layout)
    ply.plot(fig, filename="KIC" + kepID + filename + ".html")


def make_trend(x_axis: list, x_ingress_start: float, x_ingress_end: float,
               max_y: float, x_egress_start: float, x_egress_end: float,
               depth_: float, gradient_: float) -> list:
    """Creates piece-wise function for phase-folded plot"""

    y_trend = np.empty(len(x_axis))
    for i in range(0, len(x_axis)):
        if x_axis[i] < x_ingress_start:
            y_trend[i] = 0.0

        elif x_axis[i] > x_ingress_start and x_axis[i] < x_ingress_end:
            y_trend[i] = 0.0 + gradient_ * (x_axis[i] - x_ingress_start)

        elif x_axis[i] > x_ingress_end and x_axis[i] < x_egress_start:
            y_trend[i] = -depth_

        elif x_axis[i] > x_egress_start and x_axis[i] < x_egress_end:
            y_trend[i] = 0.0 - gradient_ * (x_axis[i] - x_egress_end)

        elif x_axis[i] > x_egress_end:
            y_trend[i] = 0.0

    return y_trend


def variance_squared(y_axis:list) -> float:
    """Estimates the noise within in the phase-folded data"""

    mean = sum(y_axis) / len(y_axis)
    var_squared= sum((i - mean) ** 2 for i in y_axis) / len(y_axis)
        
    return var_squared



def chi_squared_test(x_axis: list, y_axis: list, model: list,
                     variance_squared: float) -> float:
    """Finds the error estimate for your line of best-fit"""

    m = len(x_axis)
    chi_square = 0.0
    for i in range(0, m):
        chi_square = chi_square + \
                     (1.0 / m) * ((y_axis[i] - model[i]) ** 2.0 / variance_squared)

    return chi_square


def planet_parameters(t_egress_end: float, t_ingress_start: float,
                      t_egress_start: float, t_ingress_end: float,
                      period: float, _depth: float, r_star_input: float,
                      g_input: float) -> float:
    """Calculates planet parameters based on values associated with the 
    best-fit solution"""

    au = 1.496 * 1.0e+11
    m_solar = 1.989 * 1.0e30
    r_earth = 6.371 * 1.0e6
    g_constant = 6.67408 * 1.0e-11
    r_jupiter = 6.9911e7
    r_neptune = 2.4622e7
    T_dur = abs(t_egress_end - t_ingress_start) * 24.0 * 3600.0  
    T_dur_full = abs(t_egress_start - t_ingress_end) * 24.0 
    p_seconds = period * 24.0 * 3600.0
    r_planet = math.sqrt(_depth) * r_star_input
    order = r_planet / r_earth
    j_order = r_planet / r_jupiter
    n_order = r_planet / r_neptune
    m_star = ((10 ** g_input) * 0.01 * (r_star_input ** 2.0)) / g_constant
    m_ratio = m_star / m_solar
    a = ((g_constant * m_star * (p_seconds ** 2.0)) /
         (4.0 * (math.pi ** 2.0))) ** (1.0 / 3.0)
    a_au = a / au
    print(r_planet)
    print(T_dur)
    print(r_star_input)
    
    section1= (r_star_input + r_planet) ** 2.0
    section2= T_dur * math.pi / p_seconds
    section3= (section1)-(a * math.sin(section2))** 2
    b = math.sqrt(abs(section3)) / r_star_input 
    
    inc = math.acos(b * r_star_input / a) * (180.0 / math.pi)
    t_dur_hrs = T_dur / 3600.0
    r_star_solar_units = r_star_input / 7.0e8

    return r_star_solar_units, t_dur_hrs, T_dur_full, inc, b, a_au, m_ratio, \
           order, j_order, n_order


def main():
    kepID = input("Enter the KepID for your chosen light curve: ").strip()
    input("Rename your chosen light curve file to " + kepID + \
          ".tbl in the folder " + os.getcwd() + " then press Enter to proceed.")

    try:
        x, y = np.loadtxt(kepID + ".tbl",
                          usecols=(1, 2), skiprows=3, unpack=True)
        print("Plot loaded successfully.")
    except OSError:
        path = os.path.join(os.getcwd() + kepID + ".tbl")
        print("Data not found:", path)
        return

    r_star = get_float("Enter the stellar radius [Solar radii] for KIC"
                       + kepID + ": ") * 7.0e8
    g = get_float("Enter the stellar surface "
                  "gravity [log10(cm/s**2)]: ")

    repeat = True
    while repeat:

        original_data_trace = make_trace(x, y, 'blue')
        make_plot([original_data_trace], kepID, "Time (Julian days)")
        print("Plot saved as: KIC" + kepID + ".html")
    
        
        input("Zoom in on a section of the light curve which contains "
              "at least 10 transits without gaps (where possible). "
              "Press Enter to continue.")

        n, first_midtransit, last_midtransit = get_float_list(
            "Enter the number of transits in your region, "
            "the mid-transit time for the first transit and "
            "mid-transit time of the last transit, seperated "
            "by commas: ", 3)
        
        P = (abs(last_midtransit - first_midtransit)) / (n - 1.0)
        P = round(P, 4)
        print()
        print("Based on these values, the orbital period is", P, "days.")
        print("Loading phase-folded light curve.")
        y2 = y
        offset = (first_midtransit - (P / 2.0))
        t = x - offset


        t = (t % P) - (P / 2.0)
        t, y2 = zip(*sorted(zip(t, y2)))

        phase_fold_trace = make_trace(t, y2, 'blue')
        make_plot([phase_fold_trace], kepID, "Phase (days)",
                  "_phase_fold")
        print("Plot saved as: KIC" + kepID + "_phase_fold.html")
        repeat = ask_y_n("Would you like to reanalyse the light curve?")
        if repeat:
            input("Please close all plots and then press Enter to continue.")

    repeat_chi = True
    while repeat_chi:
        print()
        print("In order to find a best-fit model for your phase folded plot, "
              "the program will scan over a range of ingress start times, "
              "ingress end times and minimum fluxes. You will need to use "
              "the axes to determine the end points for these ranges and the "
              "transit time midpoint. Press enter to continue.")

        axis_of_symmetry = get_float("Enter the time of the midpoint of the "
                                     "transit: ")

        t_ingress_start_min, t_ingress_start_max = get_float_list(
            "Enter minimum and maximum time for ingress start, "
            "seperated by a comma: ", 2)

        t_ingress_end_min, t_ingress_end_max = get_float_list(
            "Enter a minimum and maximum time for ingress end, "
            "seperated by a comma: ", 2)

        min_flux_min, min_flux_max = get_float_list(
            "Enter a minimum and maximum values for the avarage "
            "minimum flux seperated by a comma: ", 2)

        t_ingress_start_array = make_array(t_ingress_start_min,
                                           t_ingress_start_max)
        t_ingress_end_array = make_array(t_ingress_end_min, t_ingress_end_max)
        min_flux_array = make_array(min_flux_min, min_flux_max)
        debug()
        debug(t_ingress_start_array)
        debug(t_ingress_end_array)
        debug(min_flux_array)
        debug()

        varsqr = variance_squared(y2)
        product = itertools.product(t_ingress_start_array,
                                    t_ingress_end_array, min_flux_array)

        best_chi_square = 1
        best_chi_inputs = []
        best_y_model = np.empty(len(t))
        best_t_ingress_start = 0.0
        best_t_ingress_end = 0.0
        best_min_flux = 0.0
        best_t_egress_start = 0.0
        best_t_egress_end = 0.0
        
        print()
        print("Finding a best-fit model based on these ranges.")
        start = time.time()
        for i, combo in enumerate(product):
            t_ingress_start = combo[0]
            t_ingress_end = combo[1]
            min_flux = combo[2]
            max_flux = 0.0

            t_egress_start = axis_of_symmetry + (axis_of_symmetry - t_ingress_end)
            t_egress_end = axis_of_symmetry + (axis_of_symmetry - t_ingress_start)

            depth = abs(max_flux - min_flux)
            gradient = (min_flux - max_flux) / (t_ingress_end -
                                                t_ingress_start)

            y_model = make_trend(t, t_ingress_start, t_ingress_end, max_flux,
                                 t_egress_start, t_egress_end, depth, gradient)
            chisqr = chi_squared_test(t, y2, y_model, varsqr)

            if i % 100 == 0:
                print(".", end="")

         
            if chisqr < best_chi_square:
                debug("A model was found with a Chi-square value of:", chisqr)
                best_chi_square = chisqr
                best_chi_inputs = combo
                best_t_ingress_start = t_ingress_start
                best_t_ingress_end = t_ingress_end
                best_min_flux = min_flux
                best_t_egress_start = t_egress_start
                best_t_egress_end = t_egress_end
                best_y_model = y_model
                

        end = time.time()
        diff = int(round(end - start))
        debug("Variance:", varsqr)
        print()
        print("Lowest Chi-Square value of", best_chi_square, "was found in",
              diff, "seconds.")
        debug("The associated ingress start time, ingress end time and "
              "minimum flux were found:", best_chi_inputs)
        debug("The associated egress start time, egress end time were found",
              best_t_egress_start, best_t_egress_end)

        best_depth = abs(max_flux - best_min_flux)

        debug("The associated depth:", depth)
        trace_2 = make_trace(t, y2, 'blue', 'Phase Folded Data')
        trace_3 = make_trace(t, best_y_model, 'red', 'Model')
        make_plot([trace_2, trace_3], kepID, "Phase (days)", "_trend")
        print("Plot saved as: KIC" + kepID + "_trend.html")
        r_star_solar_units, t_dur_hrs, T_dur_full, inc, b, a_au, m_ratio, \
        order, j_order, n_order = planet_parameters \
            (best_t_egress_end, best_t_ingress_start, best_t_egress_start,
             best_t_ingress_end, P, best_depth, r_star, g)

        print()
        print("The following parameters for the KIC" + kepID + " system"
                                                               " were found:")
        print("Planetary radius[Earth radii]: ", round(order, 2))
        print("Planetary radius[Jupiter radii]: ", round(j_order, 2))
        print("Planetary radius[Neptune radii]: ", round(n_order, 2))
        print("Period[days]", P)
        print("Stellar mass[solar mass]: ", round(m_ratio, 3))
        print("Orbit semi-major axis[au]: ", round(a_au, 4))
        print("Impact parameter: ", round(b, 3))  
        print("Inclination [deg]: ", round(inc, 2))
        print("Total Transit duration[hrs]: ", round(t_dur_hrs, 5))
        print("Full Transit duration[hrs]: ", round(T_dur_full, 5))
        repeat_chi = ask_y_n("Do you want to change the range of inputs "
                             "for your model?")

    save_array = np.empty(shape=[0, 19])

    if ask_y_n("Do you want to save the inputs and parameters?"):
        filename_ = input("Enter a filename: ") + ".csv"
        header = ["KepID"
            , "Period[days]"
            , "Planetary Radius[Earth radii]"
            , "Planetary Radius[Jupiter radii]"
            , "Planetary Radius[Neptune radii]"
            , "Stellar Mass[solar mass]"
            , "Orbit Semi-Major Axis[au]"
            , "Impact Parameter"
            , "Inclination[deg]"
            , "Total Transit duration[hrs]"
            , "Full Transit duration[hrs]"
            , "Stellar Radius[Solar radii]"
            , "Stellar Surface Gravity[log10(cm/s**2)]"
            , "Ingress Start Time"
            , "Ingress End Time"
            , "Min flux"
            , "Egress Start Time"
            , "Egress End Time"
            , " Chi-Square"
                  ]
        headers = ",".join(header)
        values = np.array([[int(kepID), P, order, j_order, n_order, m_ratio,
                            a_au, b, inc, t_dur_hrs, T_dur_full,
                            r_star_solar_units, g, best_t_ingress_start,
                            best_t_ingress_end, best_min_flux,
                            best_t_egress_start, best_t_egress_end,
                            best_chi_square]])

        save_array = np.append(save_array, np.stack(values), axis=0)

        try:
            np.savetxt(filename_, save_array, delimiter=",",
                       header=headers)
            print("Saved to: ", filename_)
        except PermissionError:
            print("Permission denied trying to save file.")
            print("Please close the file first and run the program again.")
            return

        if ask_y_n("Would you like to open " + filename_ + " now?"):

            success = os.system('start "" ' + filename_)
            if success != 0:
                print("Could not open file.")


if __name__ == "__main__":
    main()
