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
Module DetDim
! This is used for the MCLR code
! contains all PARAMETERS defining LUCIA
INTEGER, PARAMETER ::MXPIRR = 20
INTEGER, PARAMETER ::MXPOBS = 20
INTEGER, PARAMETER ::MXPR4T = 10
INTEGER, PARAMETER ::MXINKA = 200   ! Resolution of identity
INTEGER, PARAMETER ::MXPORB = 500   ! Maximum number of orbitals
INTEGER, PARAMETER ::MXPXOT = 9
INTEGER, PARAMETER ::MXPXST = 100
INTEGER, PARAMETER ::MXPSHL = 100
INTEGER, PARAMETER ::MXPL = 20
INTEGER, PARAMETER ::MXPXT = 25
INTEGER, PARAMETER ::MXPICI = 30
INTEGER, PARAMETER ::MXPSTT = 2500
INTEGER, PARAMETER ::MXPCSM = 20
INTEGER, PARAMETER ::MXPCTP = 30
INTEGER, PARAMETER ::MXCNSM = 8
!
INTEGER, PARAMETER ::MXPWRD = 2000000
INTEGER, PARAMETER ::MXNMS = 5
INTEGER, PARAMETER ::MTYP = 30

!. Note : MXPNGAS = MXPR4T+6 !!
!. Required in order to handle GAS and RAS within /LUCINP/
INTEGER, PARAMETER ::MXPNGAS = 3
INTEGER, PARAMETER::MXPNSMST = 8
!. Largest allowed division of space for perturbation operator
INTEGER, PARAMETER ::MXPPTSPC=20
End Module DetDim
