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
SubRoutine Init_Int_Options()
use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
use Int_Options, only: Thize, W2Disc, PreSch, Disc_Mx, Disc, Quad_ijkl
use Constants, only: Zero
Implicit None

!     Set variables in module Int_Options

DoIntegrals=.True.  ! Default value
DoFock=.False.      ! Default value
FckNoClmb=.False.   ! Default value
FckNoExch=.False.   ! Default value
Thize=Zero          ! Default value
W2Disc=.False.      ! Default value
PreSch=.True.       ! Default value
Disc_Mx=Zero        ! Default value
Disc=Zero           ! Default Value
Quad_ijkl=Zero      ! Default Value

End SubRoutine Init_Int_Options
