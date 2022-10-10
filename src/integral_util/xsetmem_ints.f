************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine xsetmem_ints(mem)
      use k2_arrays, only: Sew_Scr
      implicit real*8 (a-h,o-z)
#include "status.fh"
#include "stdalloc.fh"
c
      If (XMem_Status.eq.Active) Then
         Call WarningMessage(2,
     &               'External handling of scratch already active!')
         Call Abend()
      End If
      Mem_ = Mem
C Avoid using up all available memory
      Call mma_maxDBLE(MemMax)
      If (MemMax-Mem_.lt.1000) Then
        If (Mem_.gt.1000) Mem_=Mem_-1000
      End If
C     Write (6,*) 'xsetmem_ints: External allocate:',Mem_
      Call mma_allocate(Sew_Scr,Mem_,Label='Sew_Scr')
      XMem_Status=Active
C     Call mma_MaxDBLE(nu)
C     Write (6,*) 'xsetmem_ints: External allocate left to allocate:',nu
*
      Return
      End
