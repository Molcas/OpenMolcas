!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! THERMOstat = 0  nothing (default)
!              1  microcanonical ensemble (Andersen Thermostat)
!              2  canonical ensemble (Nose-Hoover Chain of Thermostats)
! VELOcities = 0  zero start velocities (default)
!              1  reads in the velocities (Bohr/a.u.) from 'velocity.xyz'
!              2  reads in the mass-weighted velocities (Bohr/a.u.) from 'velocity.xyz'
!              3  Maxwell-Boltzmann distribution at temperature T (Bohr/a.u.)
! POUT = x Number of nuclear coordinates to project out from the dynamics
!          and therefore of files 'out.00x.xyz' to read
! PIN = x Number of nuclear coordinates to keep in for the dynamics
!          and therefore of files 'in.00x.xyz' to read
!
! DT      - Time step
! Etot0   - Reference total energy for microcanonical ensemble
! Restart - Restart time of the molecular dynamics simulation
! iPrint  - (global) print level

module Dynamix_Globals

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: THERMO, VELO, POUT, PIN, iPrint
real(kind=wp) :: DT, RESTART, TEMP
logical(kind=iwp) :: lH5Restart
integer(kind=iwp), parameter :: SILENT = 0, TERSE = 1, USUAL = 2, VERBOSE = 3, DEBUG = 4, INSANE = 5
integer(kind=iwp), parameter :: nh = 6, iQ1 = 1, iQ2 = 2, iX1 = 3, iX2 = 4, iVx1 = 5, iVx2 = 6
integer(kind=iwp), parameter :: VelVer = 1, VV_First = 2, VV_Second = 3, Gromacs = 4

public :: DT, lH5Restart, iPrint, PIN, POUT, THERMO, TEMP, RESTART, VELO
public :: SILENT, TERSE, USUAL, VERBOSE, DEBUG, INSANE
public :: nh, iQ1, iQ2, iX1, iX2, iVx1, iVx2
public :: VelVer, VV_First, VV_Second, Gromacs

#ifdef _HDF5_
integer(kind=iwp) :: dyn_fileid, dyn_time, dyn_dt, dyn_etot, dyn_etot0, dyn_vel, dyn_nh, dyn_mass, dyn_geom
character(len=180) :: File_H5Res
public :: dyn_fileid, dyn_time, dyn_dt, dyn_etot, dyn_etot0, dyn_vel, dyn_nh, dyn_mass, dyn_geom, File_H5Res
#endif

end module Dynamix_Globals
