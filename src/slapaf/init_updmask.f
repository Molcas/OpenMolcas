************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Roland Lindh                                     *
************************************************************************
      Subroutine Init_UpdMask(Curvilinear, Redundant, nsAtom, nInter)
      Use NewH_mod
      Implicit None
#include "stdalloc.fh"
      Integer, Allocatable:: IsMM(:)
      Logical Curvilinear, Redundant
      Integer nsAtom, nInter
*
      Integer iAtom, i, nAtMM
*
*---- If redundant Cartesians and no symmetry, use unit matrix for MM atoms
*
      If (Redundant.and.(.not.Curvilinear).and.
     &    (3*nsAtom.eq.nInter)) Then
         Call mma_allocate(UpdMask,nInter,label="UpdMask")
         Call mma_allocate(IsMM,nsAtom,Label='IsMM')
         Call MMCount(nsAtom,nAtMM,IsMM)
         Do iAtom=1,nsAtom
            If (IsMM(iAtom).eq.1) Then
               Do i=1,3
                 UpdMask((iAtom-1)*3+i)=1
               End Do
            Else
               Do i=1,3
                 UpdMask((iAtom-1)*3+i)=0
               End Do
            End If
         End Do
         Call mma_deallocate(IsMM)
      End If
*
      Return
      End Subroutine Init_UpdMask
