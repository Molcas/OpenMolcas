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

module DetDim
! This is used for the MCLR code
! contains all PARAMETERS defining LUCIA

integer, parameter :: MXPIRR = 20
integer, parameter :: MXPOBS = 20
integer, parameter :: MXPR4T = 10
integer, parameter :: MXINKA = 200   ! Resolution of identity
integer, parameter :: MXPORB = 500   ! Maximum number of orbitals
integer, parameter :: MXPXOT = 9
integer, parameter :: MXPXST = 100
integer, parameter :: MXPSHL = 100
integer, parameter :: MXPL = 20
integer, parameter :: MXPXT = 25
integer, parameter :: MXPICI = 30
integer, parameter :: MXPSTT = 2500
integer, parameter :: MXPCSM = 20
integer, parameter :: MXPCTP = 30
integer, parameter :: MXCNSM = 8
integer, parameter :: MXPWRD = 2000000
integer, parameter :: MXNMS = 5
integer, parameter :: MTYP = 30

! Note : MXPNGAS = MXPR4T+6 !!
! Required in order to handle GAS and RAS within /LUCINP/
integer, parameter :: MXPNGAS = 3
integer, parameter :: MXPNSMST = 8
! Largest allowed division of space for perturbation operator
integer, parameter :: MXPPTSPC = 20

end module DetDim
