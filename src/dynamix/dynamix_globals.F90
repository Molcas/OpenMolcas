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
character(len=180) :: File_H5Res

integer(kind=iwp) :: dyn_fileid, dyn_time, dyn_dt, dyn_etot, dyn_etot0, dyn_vel, dyn_nh, dyn_mass, dyn_geom

integer(kind=iwp), parameter :: SILENT = 0, TERSE = 1, USUAL = 2, VERBOSE = 3, DEBUG = 4, INSANE = 5

public :: DT, File_H5Res, iPrint, lH5Restart, PIN, POUT, THERMO, TEMP, RESTART, VELO
public :: SILENT, TERSE, USUAL, VERBOSE, DEBUG, INSANE
public :: dyn_fileid, dyn_time, dyn_dt, dyn_etot, dyn_etot0, dyn_vel, dyn_nh, dyn_mass, dyn_geom

end module Dynamix_Globals
