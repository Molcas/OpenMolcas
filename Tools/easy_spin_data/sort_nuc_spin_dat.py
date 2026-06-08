#!/usr/bin/env python3
# ***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2026, Dong Q. Le                                       *
# ***********************************************************************

import os

import requests

# EasySpin isotope data stored at:
# MOLCAS_DIR/Tools/spin_data/isotopedata.txt
SCRIPT_DIR       = os.path.dirname(os.path.abspath(__file__))
TOOLS_DIR        = "/".join(SCRIPT_DIR.split("/")[:-1])
MOLCAS_DIR       = "/".join(SCRIPT_DIR.split("/")[:-2])
DATA_DIR         = MOLCAS_DIR + "/data"
MOLCAS_SPIN_DATA = DATA_DIR + "/spin_data.src"
EASYSPIN_DATA    = SCRIPT_DIR + "/isotopedata.txt"

print("------------")
print("Tools dir  : ", TOOLS_DIR)
print("Data dir   : ", DATA_DIR)
print("------------")
print()

# Raw GitHub URL for the file
GITHUB_URL = (
    "https://raw.githubusercontent.com/"
    "StollLab/EasySpin/main/easyspin/private/isotopedata.txt"
)


def ensure_isotopedata():
    """Download/refresh isotopedata.txt from GitHub."""
    # If file does not exist, or you always want the latest, download it
    if not os.path.exists(EASYSPIN_DATA):
        print("File not found, downloading...")
        download_file()
    else:
        print("File already exists, not downloading. ", EASYSPIN_DATA)
        print("If you want to update EasySpin database, simply delete isotopedata.txt file then re-run.")
        print("Please check the data carefully (format, g-factor or magnetic moment, units,...) before generating a new database.")


def download_file():
    resp = requests.get(GITHUB_URL, timeout=30)
    resp.raise_for_status()  # raise error if download failed
    with open(EASYSPIN_DATA, "wb") as f:
        f.write(resp.content)
    print(f"Saved latest file to {EASYSPIN_DATA}")


ensure_isotopedata()

isotopes = []

with open(EASYSPIN_DATA, "r") as isotopedatatxt:
    line = isotopedatatxt.readline()
    while line:
        if line.startswith("%"):
            line = isotopedatatxt.readline()
            continue
        else:
            data = []
            columns = line.split()
            Z = int(columns[0])
            A = int(columns[1])
            is_stable = columns[2]
            symbol = columns[3]
            spin = float(columns[5])
            magnetic_moment = float(columns[6])
            # Note: magnetic_moment = gfactor * spin
            # In case magnetic_moment = 0, spin is also zero
            # The assignment gfactor = magnetic_moment is made to prevent dividing by zero.
            gfactor = magnetic_moment
            if spin != 0.0:
                gfactor = gfactor / spin
            abundance = float(columns[7])
            electric_quad = columns[8]
            isotopes.append(
                [
                    Z,
                    A,
                    abundance,
                    spin,
                    gfactor,
                    is_stable,
                    electric_quad,
                    symbol,
                ]
            )

        line = isotopedatatxt.readline()

# Sort by Z ascending, abundance descending
#                                                  AtNumb : Col 1     Abundance : Col 3
isotopes = sorted(isotopes, key=lambda columns: (int(columns[0]), -float(columns[2])))

# Create element map [begin_index, end_index] in database
iElem = 0
Elements = []
iIso = 0
numb_isotopes = len(isotopes)
firstIsotope_idx = []
lastIsotope_idx = []
old_Elem_Symbol = ""
while iIso != numb_isotopes:
    iIso = iIso + 1
    if isotopes[iIso - 1][7] != old_Elem_Symbol:
        Elements.append(isotopes[iIso - 1][7])
        iElem = iElem + 1
        firstIsotope_idx.append(iIso)
        lastIsotope_idx.append(iIso - 1) if (old_Elem_Symbol != "") else None
        old_Elem_Symbol = isotopes[iIso - 1][7]

    if iIso == numb_isotopes:
        lastIsotope_idx.append(iIso)

num_elements = len(firstIsotope_idx)

print("\n\n\n")
print("PARSING SUMMARY")
print("---------------")
print("Number of Elements: ", num_elements)
print("Number of Isotopes: ", numb_isotopes)


with open(MOLCAS_SPIN_DATA, "w") as spin_dat_file:
    spin_dat_file.write("""!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Dong Q. Le                                       *
!***********************************************************************
! Data derived from EasySpin (https://easyspin.org/)                   *
! Copyright (c) 2026                                                   *
! Licensed under the MIT License.                                      *
! Please see the file readme_easyspin_license.md in the                *
! Tools/spin_data directory for licence details.                       *
!***********************************************************************


! HEADER FORMAT---------------------------------------------------------
!
!>>> Element ___ Z: ___ Isotope Record: ___ - ___ TOTAL: ___
!           13-15  20-22               40-42 46-48     57-59
!
!Rec   AtNumb   MassNumb      Abundance   NucSpin    Nuc-Gfactor   Stable
!1-4     8-13      17-24          28-39     43-49          54-64    71-73
!
! Note: stable -, radioactive *
!-----------------------------------------------------------HEADER FORMAT
\n\n\n
""")

    spin_dat_file.write(
        f"@ NUMBER_OF_ELEMENTS {num_elements:>3d}             NUMBER_OF_ISOTOPES {numb_isotopes:>3d}\n\n"
    )

    for iElem in range(num_elements):
        FirstRec = firstIsotope_idx[iElem]
        LastRec = lastIsotope_idx[iElem]
        Element = Elements[iElem]
        TotalRec = LastRec - FirstRec + 1
        spin_dat_file.write(f"""
>>> Element {Element:<3}   Z: {iElem + 1:>3d}       Records: {FirstRec:>3d} - {LastRec:>3d}    TOTAL ISOTOPES:  {TotalRec:>3d}
!Rec   AtNumb   MassNumb      Abundance   NucSpin    Nuc-Gfactor   Stable
""")

        for iRec in range(FirstRec - 1, LastRec):
            AtNumb = isotopes[iRec][0]
            MassNumb = isotopes[iRec][1]
            Abundance = isotopes[iRec][2]
            NucSpin = isotopes[iRec][3]
            GFactor = isotopes[iRec][4]
            if isotopes[iRec][5] == "-":
                Stable = ".T."
            else:
                Stable = ".F."
            spin_dat_file.write(
                f""" {iRec + 1:>3d}   {AtNumb:>6d}   {MassNumb:>8d}   {Abundance:>12.8f}   {NucSpin:>7.1f}   {GFactor:>12.8f}   {Stable:>6}\n"""
            )

print("This tool has created the file spin_dat.src at:")
print(MOLCAS_SPIN_DATA)
print()
print()
print(
    "NOTE: [OpenMOLCAS] Default hyperfine spin database includes an additional record for 141Ce. REF: 10.61092/iaea.zhrz-3g5j"
)
print(
    "      If you download a new isotope data from EasySpin, 141Ce might not be included."
)
