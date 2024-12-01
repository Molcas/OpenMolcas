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
      Subroutine RdDens(d1,nd1,d2,nd2)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One, Two
      use files_mclr, only: LuJob
      Implicit None
      Integer nd1, nd2
      Real*8 D1(nd1),d2(nd2)
#include "Input.fh"
#include "SysDef.fh"
      Real*8 rdum(1)
      Real*8, Allocatable:: G2tt(:), D2t(:), D1t(:)
      Integer jDisk,i,j,iB,jB,iDij,kB,lB,iDkl,iIJKL
      Real*8 W,Fact
*                                                                      *
************************************************************************
*                                                                      *
      Integer itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      d1(:)=Zero
      Call mma_allocate(G2tt,nd2,Label='G2tt')
      Call mma_allocate(D2t,nd2,Label='D2t')
      Call mma_allocate(D1t,nd1,Label='D1t')
      G2tt(:)=Zero
      jDisk = ITOC(3)
      Do i=1,lroots
         W=Zero
         Do j=1,nroots
            If (iroot(j).eq.i) W=Weight(j)
         End Do
         Call dDaFile(LUJOB ,2,D1t,nd1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,nd1,jDisk)
         Call dDaFile(LUJOB ,2,D2t,ND2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ND2,jDisk)
         If (W.ne.Zero) Then
            call daxpy_(nd2,w,D2t,1,G2tt,1)
            call daxpy_(nd1,w,D1t,1,D1,1)
         End If
!         Write(*,*) i,w,"LUJOB",LUJOB
!        Call Triprt('D',' ',D1,ntash)
      End Do
      Call Put_dArray('D2av',G2tt,nd2)
      Call Put_dArray('D1av',D1,nd1)
*
      Do iB=1,ntash
       Do jB=1,iB
        iDij=iTri(ib,jB)
        Do kB=1,ib
         Do lB=1,kB
          iDkl=iTri(kB,lB)
          fact=One
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
          iijkl=itri(iDij,iDkl)
          D2(iijkl)=Fact*G2tt(iijkl)
         End Do
        End Do
       End Do
      End Do
c
      Call mma_deallocate(G2tt)
      Call mma_deallocate(D2t)
      Call mma_deallocate(D1t)
      End Subroutine RdDens
