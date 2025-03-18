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
      Subroutine mkp1(nEX,lst,rMat,rdiag)
      use Constants, only: Zero
      use negpre, only: LuCIV, P1
      use stdalloc, only: mma_allocate, mma_deallocate
      use input_mclr, only: lRoots,nConf,ERASSCF
      Implicit None
      Integer nEX
      Integer lst(nex)
      Real*8 rMat(*),rdiag(*)
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)

      Integer i,j,k,l,itri,idisk,jDisk,kk,ll
      Real*8 rtmp

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Call mma_allocate(TMP1,nconf,Label='Tmp1')
      Call mma_allocate(TMP2,nconf,Label='Tmp2')

      idisk=0
      Do i=1,lroots
       Call dDaFile(LuCIV,2,Tmp1,nconf,iDisk)
       jdisk=0
       Do j=1,i
         Call dDafile(luciv,2,Tmp2,nconf,jDisk)
         rTmp=Zero
         Do k=1,nex
          do l=1,nex
           kk=lst(k)
           ll=lst(l)
           rtmp=rtmp+Tmp1(kk)*Tmp2(ll)*rmat(itri(k,l))
          End Do
         End Do
         Do k=1,nconf
          rtmp=rtmp+Tmp1(k)*Tmp2(k)*rdiag(k)
         End Do
         If (i.eq.j) rtmp=rtmp-ERASSCF(1)
         Do  k=1,nEx
          kk=lst(k)
          rtmp=rtmp-Tmp1(kk)*Tmp2(kk)*(rdiag(kk+1)-ERASSCF(1))
         End Do
         P1(itri(i,j))=rtmp
        End Do
       End Do

       Call mma_deallocate(TMP2)
       Call mma_deallocate(TMP1)

       end Subroutine mkp1
