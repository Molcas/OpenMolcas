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
      Subroutine Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                 jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                 iCmpi,jCmpj,iShell,jShell,
     &                 mDCRij,ipDij,ipDDij,mDij,nIrrep,DeDe,nDeDe)
      Implicit Real*8 (a-h,o-z)
      Real*8 DeDe(nDeDe)
#include "WrkSpc.fh"
*
      If (nIrrep.eq.1) Then
         ii1=0
         ii2=1
         ii3=0
         jj1=0
         jj2=1
         jj3=0
      Else
         ii1=iBasi
         ii2=iBasAO
         ii3=iBasn
         jj1=jBasj
         jj2=jBasAO
         jj3=jBasn
      End If
      If (mDCRij.ne.0) Then
         If (iShell.ge.jShell) Then
            i1=ii1
            i2=ii2
            i3=ii3
            j1=jj1
            j2=jj2
            j3=jj3
         Else
            i1=jj1
            i2=jj2
            i3=jj3
            j1=ii1
            j2=ii2
            j3=ii3
         End If
         If (iBasi.eq.iBsInc.and.jBasj.eq.jBsInc) Then
            ipDDij=ipDij
         Else
            Call Picky(DeDe(ipDij),i1,j1,iPrimi*jPrimj,
     &                 iCmpi*jCmpj,mDCRij,
     &                 i2,i2+i3-1,j2,j2+j3-1,DeDe(ipDDij))
         End If
      End If
      mDij = (ii3*jj3+1)*iCmpi*jCmpj + iPrimi*jPrimj + 1
*
      Return
      End
