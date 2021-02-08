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
      use negpre
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "stdalloc.fh"
#include "incdia.fh"
      Real*8 rdia(*)
      Real*8, Allocatable:: TMP1(:), TMP2(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Call mma_allocate(TMP1,nconf,Label='TMP1')
      Call mma_allocate(TMP2,nconf,Label='TMP2')

      idisk=0
      Do i=1,lroots
       jdisk=idisk
       Call dDaFile(LuCIV,2,Tmp1,nconf,iDisk)
       Call ExpHinvv(rdia,Tmp1,Tmp1,0.0d0,1.0d0)
       Do j=i,lroots
         Call dDafile(luciv,2,Tmp2,nconf,jDisk)
         p1INV(itri(i,j))=DDOT_(nconf,Tmp2,1,Tmp1,1)
       End Do
      End Do

      Call mma_deallocate(TMP1)
      Call mma_deallocate(TMP2)

      Return
      end
