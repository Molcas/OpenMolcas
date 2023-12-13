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
Module Disp
Private
#include "Molcas.fh"
Logical TRSymm, lEq, Dirct(MxAtom*3), l2DI, HF_Force

Character(LEN=LENIN6) ChDisp(MxAtom*3)

Integer IndDsp(MxAtom,0:7), iSkal(MxBas), InxDsp(MxAtom,3),&
        lDisp(0:7), IndxEq(MxAtom*3), ipAM, nTR,           &
        mult_Disp(MxAtom*3)

Real*8    CutGrd, Disp_Fac(3,0:7,MxAtom)

Public :: TRSymm, lEq, Dirct, l2DI, HF_Force, ChDisp, IndDsp, iSkal, InxDsp,&
          lDisp, IndxEq, ipAM, nTR, mult_Disp,  CutGrd, Disp_Fac
End Module Disp
