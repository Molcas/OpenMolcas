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
      Subroutine mkcipre
      use negpre
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "Pointers.fh"
#include "Input.fh"

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      irec(i,j)=i+(j-1)*2*lroots


      Call mma_allocate(SS,4*lroots**2,Label='SS')
      DO I=1,lroots
       DO J=1,lroots
        SS(irec(2*i-1,2*j-1))=P1(itri(i,j))
       End Do
      End Do
      DO I=1,lroots
        SS(irec(2*i-1,2*i-1))= SS(irec(2*i-1,2*i-1))+ERAS(I)-ERASSCF(1)
        SS(irec(2*i,2*i-1))=-1.0d0
        SS(irec(2*i-1,2*i))=-1.0d0
      End Do
      SS(irec(2*lroots-1,2*lroots-1))=
     &     SS(irec(2*lroots-1,2*lroots-1))+1.0d0
      Call MatInvert(SS,2*lroots)
      DO I=1,lroots
       DO J=1,lroots
        SS(irec(2*i-1,2*j-1))= SS(irec(2*i-1,2*j-1))+P1INV(itri(i,j))
        SS(irec(2*i,2*j))= SS(irec(2*i,2*j))+P1(itri(i,j))
       End Do
      End Do
      DO I=1,lroots
          SS(irec(2*i,2*i-1)) = SS(irec(2*i,2*i-1)) + 1.0d0
          SS(irec(2*i-1,2*i)) = SS(irec(2*i-1,2*i)) + 1.0d0
      End Do
      Call MatInvert(SS,2*lroots)
      Call DSCAL_(4*lroots**2,-1.0d0,SS,1)
      SS(irec(2*lroots,2*lroots))= SS(irec(2*lroots,2*lroots))-1.0d0

      Return
      End
