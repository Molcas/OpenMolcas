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

module Disp

use Molcas, only: LenIn6, MxAtom
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: IndDsp(MxAtom,0:7), IndxEq(MxAtom*3), InxDsp(MxAtom,3), lDisp(0:7), mult_Disp(MxAtom*3), nTR
real(kind=wp) :: CutGrd, Disp_Fac(3,0:7,MxAtom)
logical(kind=iwp) :: Dirct(MxAtom*3), HF_Force, l2DI, lEq, TRSymm
character(len=LenIn6) :: ChDisp(MxAtom*3)

public :: ChDisp, CutGrd, Dirct, Disp_Fac, HF_Force, IndDsp, IndxEq, InxDsp, l2DI, lDisp, lEq, mult_Disp, nTR, TRSymm

end module Disp
