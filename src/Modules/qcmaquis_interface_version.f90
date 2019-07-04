!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_version

! stefan: DMRG interface version variables

  implicit none

! QCMaquis-interface version
/* #undef QCMaquis_interface_VERSION */
!
#define QCMAQUIS_INTERFACE_VERSION_MAJOR 2
/* #undef QCMAQUIS_INTERFACE_VERSION_MINOR */
#define QCMAQUIS_INTERFACE_VERSION_BUILD N/A (N/A)
#define QCMAQUIS_INTERFACE_GIT_VERSION  "N/A (N/A)"
!
! QCMaquis-interface version (full string)
#define QCMAQUIS_INTERFACE_VERSION_STRING "QCMaquis driver - version: 2.0"

character(len=200), public :: qcmaquis_interface_v = QCMAQUIS_INTERFACE_VERSION_STRING
character(len=200), public :: qcmaquis_interface_g = QCMAQUIS_INTERFACE_GIT_VERSION

end module qcmaquis_interface_version
