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
      Subroutine mkp1inv(rdia)
      use Constants, only: Zero,One
      use negpre, only: LuCIV, P1Inv
      use stdalloc, only: mma_allocate, mma_deallocate
      use input_mclr, only: lRoots,nConf
      Implicit None
      Real*8 rdia(*)
      Real*8, Allocatable:: TMP1(:), TMP2(:)
      Integer i,j,itri,iDisk,jDisk
      Real*8, External:: DDot_

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Call mma_allocate(TMP1,nconf,Label='TMP1')
      Call mma_allocate(TMP2,nconf,Label='TMP2')

      idisk=0
      Do i=1,lroots
       jdisk=idisk
       Call dDaFile(LuCIV,2,Tmp1,nconf,iDisk)
       Call ExpHinvv(rdia,Tmp1,Tmp1,Zero,One)
       Do j=i,lroots
         Call dDafile(luciv,2,Tmp2,nconf,jDisk)
         p1INV(itri(i,j))=DDOT_(nconf,Tmp2,1,Tmp1,1)
       End Do
      End Do

      Call mma_deallocate(TMP1)
      Call mma_deallocate(TMP2)

      end Subroutine mkp1inv
