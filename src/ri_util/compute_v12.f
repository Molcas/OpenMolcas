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
      SubRoutine Compute_V12(V,V12,nDim)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Real*8 V(nDim,nDim), V12(nDim,nDim)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
      Call Allocate_Work(ipVec,nDim**2)
      Call Allocate_Work(ipVTri,nDim*(nDim+1)/2)
*
      Call Compute_V12_(V,V12,Work(ipVTri),Work(ipVec),nDim)
*
      Call Free_Work(ipVTri)
      Call Free_Work(ipVec)
*
      Return
      End
      SubRoutine Compute_V12_(V,V12,VTri,Vec,nDim)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 V(nDim,nDim), V12(nDim,nDim), VTri(nDim*(nDim+1)/2),
     &       Vec(nDim,nDim)
*
      Call FZero(Vec,nDim**2)
      call dcopy_(nDim,[One],0,Vec,nDim+1)
*
      Do i = 1, nDim
         Do j = 1, i
            VTri(i*(i-1)/2+j) = V(i,j)
         End Do
      End Do
*
      Call JACOB(VTri,Vec,nDim,nDim)
*
      Call FZero(V12,nDim**2)
      Do i = 1, nDim
#ifdef _DEBUGPRINT_
         tmp=VTri(i*(i+1)/2)
         Write (6,*) 'i,tmp=',i,tmp
#endif
         tmp=Sqrt(VTri(i*(i+1)/2))
         If (tmp.lt.1.0D-90) Then
            V12(i,i)=1.0D90
         Else
            V12(i,i)=One/tmp
         End If
      End Do
*
      Call DGEMM_('N','T',
     &            nDim,nDim,nDim,
     &            1.0d0,V12,nDim,
     &            Vec,nDim,
     &            0.0d0,V,nDim)
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,Vec,nDim,
     &            V,nDim,
     &            0.0d0,V12,nDim)
*
      Return
      End
