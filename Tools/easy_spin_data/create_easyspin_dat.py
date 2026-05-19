#!/usr/bin/env python3
#***********************************************************************
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
#***********************************************************************

import os
import requests

# EasySpin isotope data stored at:
# MOLCAS_ROOT/Tools/create_easy_spin_data/isotopedata.txt
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOOLS_DIR = '/'.join(SCRIPT_DIR.split('/')[:-1])
SRC_DIR = '/'.join(SCRIPT_DIR.split('/')[:-2]) + '/src'
HFC_DATA_FILE=SRC_DIR + '/property_util/hfc_ezspin_dat.F90'
print("------------")
print("Tools dir  : ", TOOLS_DIR )
print("src dir    : ", SRC_DIR )
print("------------")
print()
LOCAL_PATH = SCRIPT_DIR + "/isotopedata.txt"


# Raw GitHub URL for the file
GITHUB_URL = (
    "https://raw.githubusercontent.com/"
    "StollLab/EasySpin/main/easyspin/private/isotopedata.txt"
)
print("")
print("")
print("")
def ensure_isotopedata():
    """Download/refresh isotopedata.txt from GitHub."""
    # If file does not exist, or you always want the latest, download it
    if not os.path.exists(LOCAL_PATH):
        print("File not found, downloading...")
        download_file()
    else:
        print("File already exists, not downloading. ",LOCAL_PATH)
        print("")
        print("If you want to update EasySpin database, simply delete isotopedata.txt file then re-run.")
        print("")
        print("")

def download_file():
    resp = requests.get(GITHUB_URL, timeout=30)
    resp.raise_for_status()  # raise error if download failed
    with open(LOCAL_PATH, "wb") as f:
        f.write(resp.content)
    print(f"Saved latest file to {LOCAL_PATH}")

ensure_isotopedata()

raw_data =[]
isotope_data= []
element_map = []

def format_fortran_float(value):
    if value == 'NaN':
        return 'ieee_value(0.0d0, ieee_signaling_nan)'
    else:
        return f"{float(value):11.8f}d0"

def format_fortran_spin(value):
    if value == 'NaN':
        return 'ieee_value(0.0d0, ieee_signaling_nan)'
    else:
        return f"{value}d0"

def format_boolean(value):
    if value == '*':
        return '.false.'
    else:
        return '.true.'


print("This tool has created the file hfc_ezspin_dat.F90 at:")
print(HFC_DATA_FILE)
print()
print()
SRC_DIR = '/'.join(os.getcwd().split('/')[:-2]) + '/src'

with open(LOCAL_PATH, "r") as data_file:
    line = data_file.readline()
    while(line):
        if line.startswith('%'):
            line = data_file.readline()
            continue
        else:
            data=[]
            parts = line.split()
            Z = parts[0]
            A = int(parts[1])
            is_stable = format_boolean(parts[2])
            symbol = parts[3]
            spin = format_fortran_spin(parts[5])
            g_factor  = format_fortran_float(parts[6])
            abundance = parts[7]
            electric_quad = format_fortran_float(parts[8])
            raw_data.append([Z, A, is_stable, spin, g_factor, abundance, electric_quad,symbol])

        line = data_file.readline()

# Sort by Z ascending, abundance descending
raw_data =sorted (
    raw_data,
    key=lambda r: (int(r[0]), -float(r[5]))
)
# Reformat data for Fortran output [Double precision values]
for row in raw_data:
    isotope_data.append([int(row[0]), row[1], row[2], row[3], row[4], format_fortran_float(row[5]), row[6], row[7]])
# Create element map [begin_index, end_index] in database
Z=1
i=1
numb_records=len(isotope_data)
elem_begin=[]
elem_begin.append([1, 1,'H'])
elem_end=[]
while (i!=numb_records):
    if (Z==isotope_data[i-1][0]):
        i+=1
    else:
        Z=isotope_data[i-1][0]
        elem_begin.append([Z, i,isotope_data[i-1][-1]])
        elem_end.append([Z-1, i-1,isotope_data[i-2][-1]])
        i+=1
elem_end.append([isotope_data[numb_records-1][0], numb_records,isotope_data[numb_records-1][-1]])
numb_elements=len(elem_begin)

with open(HFC_DATA_FILE, "w") as fortran_file:
    fortran_file.write(f'''!***********************************************************************
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
! Tools/create_easy_spin_data directory for licence details.           *
!***********************************************************************


module hfc_data
  use, intrinsic     :: ieee_arithmetic
  use stdalloc,      only: mma_allocate, mma_deallocate
  use Definitions,   only: iwp, wp
  implicit none

  integer, parameter :: numb_records = {len(isotope_data)}

  type :: easyspin_data

    integer(kind=iwp)             :: AtNumb
    integer(kind=iwp)             :: MassNum
    real(kind=wp)                 :: abundance
    real(kind=wp)                 :: nucspin
    real(kind=wp)                 :: gfactor
''')

# PLEASE UNCOMMENT LINES BELOW IF NEEDED

    # fortran_file.write(f'''    character(len=6)              :: notation
    # logical                       :: is_stable
    # character(len=4)              :: symbol
    # real(kind=wp)                 :: quadrupole''')

# PLEASE UNCOMMENT LINES ABOVE IF NEEDED

    fortran_file.write(f'''  end type easyspin_data

  type :: element_map
    integer(kind=iwp)             :: first_rec
    integer(kind=iwp)             :: last_rec
  end type element_map

  type(element_map),allocatable   :: elem_map_db(:)
  type(easyspin_data),allocatable :: ezspin_db(:)

    public :: init_isotope_data, gfac_spin_by_mass, get_first_nonzero_gfactor, &
              gfactor_by_nucspin, nucspin_by_gfactor, free_isotope_data

  contains

!======================================================================
  subroutine init_isotope_data()
    integer :: istat
    allocate(ezspin_db({numb_records}), stat=istat)

    if (istat > 0) then
      write(6,*) ' Error in allocating EasySpin ezspin_db.'
    end if

    allocate(elem_map_db({numb_elements}), stat=istat)

    if (istat > 0) then
      write(6,*) ' Error in allocating EasySpin elem_map_db.'
    end if

''')
    for iElem in range(1,numb_elements+1):
        if iElem==elem_begin[iElem-1][0]:
          AtNumb=iElem-1
          fortran_file.write(f'''
!--> Element {elem_begin[AtNumb][2]}
      elem_map_db({iElem})%first_rec = {elem_begin[AtNumb][1]}
      elem_map_db({iElem})%last_rec  = {elem_end[AtNumb][1]}
!----------------------------''')


    for i, isotope in enumerate(isotope_data):
        fortran_file.write(f'''
!--> Record {i+1:03d}    = {isotope_data[i][1]}{isotope_data[i][7]}
    ezspin_db({i+1})%AtNumb      = {isotope_data[i][0]}
    ezspin_db({i+1})%MassNum     = {isotope_data[i][1]}
    ezspin_db({i+1})%abundance    = {isotope_data[i][5]}
    ezspin_db({i+1})%nucspin      = {isotope_data[i][3]}
    ezspin_db({i+1})%gfactor      = {isotope_data[i][4]}
''')

  # PLEASE UNCOMMENT LINES BELOW IF NEEDED
    #     fortran_file.write(f'''
    # ezspin_db({i+1})%notation     = '{isotope_data[i][1]}{isotope_data[i][7]:<6}'
    # ezspin_db({i+1})%is_stable    = {isotope_data[i][2]}
    # ezspin_db({i+1})%quadrupole   = {isotope_data[i][6]}'''
# PLEASE UNCOMMENT LINES ABOVE IF NEEDED

        fortran_file.write(f'''!----------------------------''')

    fortran_file.write(f'''

  end subroutine init_isotope_data''')

    real_subroutine = [ 'nucspin', 'gfactor','abundance','quadrupole']


    fortran_file.write(f'''
!======================================================================
  function gfactor_by_nucspin(AtNumb,nucspin) result(gfactor)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: nucspin
    real(kind=wp)                :: gfactor
    integer(kind=iwp)            :: iRec, first_rec, last_rec

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    gfactor = ieee_value(0.0d0, ieee_signaling_nan)  ! default value if not found
    do iRec=first_rec,last_rec
      if (abs(ezspin_db(iRec)%nucspin - nucspin) < 1e-5_wp) then
        gfactor = ezspin_db(iRec)%gfactor
        exit
      endif
    enddo

    if (ieee_is_nan(gfactor)) then
      gfactor = 0.0d0
      call warning_not_found(1)
    endif

    call isot_info(AtNumb, gfactor, nucspin)

  end function gfactor_by_nucspin
!======================================================================


!======================================================================
  function nucspin_by_gfactor(AtNumb,gfactor) result(nucspin)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: gfactor
    real(kind=wp)                :: nucspin
    integer(kind=iwp)            :: iRec, first_rec, last_rec

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    nucspin = ieee_value(0.0d0, ieee_signaling_nan)  ! default value if not found
    do iRec=first_rec,last_rec
      if (abs(ezspin_db(iRec)%gfactor - gfactor) < 1e-3_wp) then
        nucspin = ezspin_db(iRec)%nucspin
        exit
      endif
    enddo

    if (ieee_is_nan(nucspin)) then
      nucspin = 0.5d0
      call warning_not_found(0)
    endif

    call isot_info(AtNumb, gfactor, nucspin)

  end function nucspin_by_gfactor
!======================================================================


!======================================================================
  subroutine gfac_spin_by_mass(AtNumb,MassNum,gfactor,nucspin)
    integer(kind=iwp),intent(in) :: AtNumb, MassNum
    real(kind=wp), intent(out)   :: nucspin, gfactor
    integer(kind=iwp)            :: iRec, first_rec, last_rec

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    nucspin = ieee_value(0.0d0, ieee_signaling_nan)  ! default value if not found
    gfactor = ieee_value(0.0d0, ieee_signaling_nan)  ! default value if not found
    do iRec=first_rec,last_rec
      if (ezspin_db(iRec)%MassNum == MassNum) then
        nucspin = ezspin_db(iRec)%nucspin
        gfactor = ezspin_db(iRec)%gfactor
        exit
      endif
    enddo


    if (ieee_is_nan(nucspin)) then
      nucspin = 0.5d0
      call warning_not_found(0)
    endif

    if (ieee_is_nan(gfactor)) then
      gfactor = 0.0d0
      call warning_not_found(1)
    endif

    call isot_info(AtNumb, gfactor, nucspin, MassNum)

  end subroutine gfac_spin_by_mass
!======================================================================


!======================================================================
  subroutine get_first_nonzero_gfactor(AtNumb, gfactor,nucspin)
    integer(kind=iwp),intent(in)           :: AtNumb
    real(kind=wp), intent(out)             :: gfactor, nucspin
    integer(kind=iwp)                      :: iRec, first_rec, last_rec

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    gfactor = ieee_value(0.0d0, ieee_signaling_nan)  ! default value if not found
    do iRec=first_rec,last_rec
      if (ezspin_db(iRec)%gfactor /= 0.0_wp) then
        gfactor = ezspin_db(iRec)%gfactor
        nucspin = ezspin_db(iRec)%nucspin
        exit
      endif
    enddo

    if (ieee_is_nan(gfactor)) then
      gfactor = 0.0d0
      call warning_not_found(1)
    endif

    call isot_info(AtNumb, gfactor, nucspin)
  end subroutine get_first_nonzero_gfactor
!======================================================================


!======================================================================
  subroutine isot_info(AtNumb, gfactor, nucspin,MassNum)
    integer(kind=iwp), intent(in) :: AtNumb
    integer(kind=iwp), intent(in), optional :: MassNum
    real(kind=wp), intent(in)     :: gfactor, nucspin

    write(6,'(11X,A18)')    repeat('_', 18)
    write(6,'(11X,A18,I0)')    'INFO          Z = ', AtNumb
    if (present(MassNum)) then
      write(6,'(11X,A18,I0)')    '              A = ', MassNum
    endif
    write(6,'(11X,A18,F12.8)') '       g-factor = ', gfactor
    write(6,'(11X,A18,F5.2)')  '        nucspin = ', nucspin
    write(6,*) ''
  end subroutine isot_info
!======================================================================


!======================================================================
  subroutine warning_not_found(option)
    integer(kind=iwp), intent(in) :: option

    if(option == 0) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to nucspin = 0.5. This value might be fictitious. '
    else if(option == 1) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to g-factor = 0.0. This value might be fictitious.'
    endif
  end subroutine warning_not_found
!======================================================================


!======================================================================
  subroutine free_isotope_data()
    integer :: istat

    deallocate(ezspin_db, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin ezspin_db.'
    end if

    deallocate(elem_map_db, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin elem_map_db.'
    end if

  end subroutine free_isotope_data
''')

    fortran_file.write(f"\nend module hfc_data\n")