#
# Gabriella Masiha
#
# Generates and saves a random selection of n kepIDs with their corresponding 
# stellar surface gravity and stellar radius values. 
#
# Data should be downloaded from the NASA Exoplanet Archive: 
# https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
# 
# KepID must be the first column and the following columns must remain: 
# KepID, Stellar Surface Gravity and Stellar Radius
#
# CSV must be in the same working directory as this script
#
# Tested on Python 3.7.9 with NumPy version 1.19.2
#

import os

import numpy as np


def get_int(prompt: str):
    while True:
        try:
            ints = int(input(prompt))

        except ValueError:
            print("Invalid input. Please try again.")
            continue

        if ints <= 0:
            print("Integer must be greater than 0. Please try again.")
            continue

        return ints


def main():
    input("Rename your selection of Kepler data to 'kepler_data.csv' in the "
          "folder " + os.getcwd() + ", then press Enter to proceed.")

    rows = get_int('Please enter the number of rows starting'
                   ' with "#" and "KepID": ')
    grav_col = get_int('Please enter the column number of stellar '
                       'surface gravity (columns start at 0): ')
    rad_col = get_int('Please enter the column number of stellar '
                      'radius (columns start at 0): ')

    try:
        data = np.loadtxt(fname="kepler_data.csv",
                          delimiter=',', comments='#', skiprows=rows,
                          usecols=(0, grav_col, rad_col))
        print("Data loaded successfully.")
    except OSError:
        path = os.path.join(os.getcwd(), "kepler_data.csv")
        print("Data not found:", path)
        return

    n = get_int("Insert the number of random choices to be generated: ")

    rows = data.shape[0]
    print("Found", rows, "row(s) in Kepler data")

    random_idx = np.random.choice(rows, size=n, replace=False)
    random_selection = data[random_idx, :]

    filename = input("Selection complete. Enter a file name for your "
                     "random selection to be saved as: ") + ".csv"

    headers = "KepID,Stellar Surface Gravity[log10(cm/s**2)],Stellar Radius" \
              "[Solar radii]"
    try:
        np.savetxt(filename, random_selection, delimiter=",",
                   header=headers)
        print("Saved", n, "random choices to", filename)
    except PermissionError:
        print("Permission denied trying to save file.")
        print("Please close the file first and run the program again.")
        return

    prompt_user = input("Would you like to open " + filename + " now?"
                                                               " [y/n] ").lower()
    if prompt_user == "y":

        success = os.system('start "" ' + filename)
        if success != 0:
            print("Could not open file.")


if __name__ == "__main__":
    main()
