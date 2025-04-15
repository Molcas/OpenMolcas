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

subroutine Init_Data()
!***********************************************************************
!                                                                      *
!     Initialize variables in common blocks                            *
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MCLR_Data, only: FNCSF2SD, FnHlf2, FNHlf3, FnJob, FnMck, FnMOTRA, FnPT2, FnQDAT, FnTemp, FNTrI1, FNTrI2, FNTRi3, FNTRI4, &
                     FNTRI5, FnTwo, LUCSF2SD, LUHlf2, LUHlf3, LuJob, LuMck, LuMOTRA, LuPt2, LuQDAT, LuTemp, LUTrI1, LUTrI2, &
                     LUTRi3, LUTRI4, LUTRI5, LuTwo
use input_mclr, only: nsRot, nUserPT, State_Sym, UserP, UserT
use Constants, only: Zero, One

implicit none

!----------------------------------------------------------------------*
!     Define files (file names and unit numbers)                       *
!----------------------------------------------------------------------*
FnJob = 'JOBIPH'
LuJob = 10
FnMck = 'MCKINT'
LuMck = 15
FnPT2 = '.RLXPT2'
LuPT2 = 17
FnTemp = 'RESP'
LuTemp = 19
FnCSF2SD = 'REORD'
LuCSF2SD = 20

FnTwo = 'ORDINT'
LuTwo = 40
FnHLF2 = 'TEMP06'
LuHLF2 = 50
FnHLF3 = 'TEMP07'
LuHLF3 = 60

FNTRI1 = 'TEMP01'
LUTRI1 = 24
FNTRI2 = 'TEMP02'
LUTRI2 = 25
FNTRI3 = 'TEMP03'
LUTRI3 = 26
FNTRI4 = 'TEMP04'
LUTRI4 = 27
FNTRI5 = 'TEMP05'
LUTRI5 = 28
! files used only within the TwoStep Run of MCLR
FnQDAT = 'QDAT' ! some temporary data is stored here
LuQDAT = 29
! this file is exactly the same as TEMP01, but is not deleted
FnMOTRA = 'MOTRA'
LuMOTRA = 30

State_Sym = 1
nUserPT = 0
UserP = One
UserT(:) = Zero
nsRot = 0

!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*

end subroutine Init_Data
