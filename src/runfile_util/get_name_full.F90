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
      Subroutine Get_Name_Full(Element)
      Implicit None
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Character*2 Element(*)
      Integer nAtom, nAtMM, i
      Logical Found
      Character(len=LENIN), Allocatable :: LabMM(:)
!
      Call Get_nAtoms_All(nAtom)
      Call Get_Name_All(Element)
!
      Call Qpg_cArray('MMO Labels',Found,nAtMM)
      If (Found) Then
        nAtMM=nAtMM/LENIN
        Call mma_allocate(LabMM,nAtMM,label='MMO Labels')
        Call Get_cArray('MMO Labels',LabMM,LENIN*nAtMM)
        Do i=1,nAtMM
          Element(nAtom+i)=LabMM(i)(1:2)
          If (Element(nAtom+i)(2:2).eq.'_') Element(nAtom+i)(2:2)=' '
        End Do
        Call mma_deallocate(LabMM)
      End If
!
      Return
      End
