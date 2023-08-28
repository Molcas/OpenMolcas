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
Module Int_Options
use Definitions, only: wp, iwp
use Constants, only: Zero
Implicit None
Private
Logical(kind=iwp) :: DoIntegrals=.True.
Logical(kind=iwp) :: DoFock=.False.
Logical(kind=iwp) :: FckNoClmb=.False.
Logical(kind=iwp) :: FckNoExch=.False.
Real(kind=wp) :: ExFac=Zero
Real(kind=wp) :: Thize=Zero
Real(kind=wp) :: Disc_Mx=Zero
Real(kind=wp) :: Disc=Zero
Logical(kind=iwp) :: W2Disc=.False.
Logical(kind=iwp) :: PreSch=.True.
Integer(kind=iwp) :: iTOffs(8**3)
Real(kind=wp) :: Quad_ijkl=Zero
Integer(kind=iwp) :: Map4(4)

Public :: DoIntegrals, DoFock, FckNoClmb, FckNoExch, ExFac, Thize, W2Disc, &
          PreSch, Disc_Mx, Disc, iTOffs, Quad_ijkl, Map4

End Module Int_Options
