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
      Subroutine mkcipre()
      use Constants, only: One
      use negpre, only: SS, ERAS, P1, P1Inv
      use stdalloc, only: mma_allocate
      use input_mclr, only: lRoots,ERASSCF
      Implicit None
      Integer i,j,itri,iRec

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
        SS(irec(2*i,2*i-1))=-One
        SS(irec(2*i-1,2*i))=-One
      End Do
      SS(irec(2*lroots-1,2*lroots-1))=
     &     SS(irec(2*lroots-1,2*lroots-1))+One
      Call MatInvert(SS,2*lroots)
      DO I=1,lroots
       DO J=1,lroots
        SS(irec(2*i-1,2*j-1))= SS(irec(2*i-1,2*j-1))+P1INV(itri(i,j))
        SS(irec(2*i,2*j))= SS(irec(2*i,2*j))+P1(itri(i,j))
       End Do
      End Do
      DO I=1,lroots
          SS(irec(2*i,2*i-1)) = SS(irec(2*i,2*i-1)) + One
          SS(irec(2*i-1,2*i)) = SS(irec(2*i-1,2*i)) + One
      End Do
      Call MatInvert(SS,2*lroots)
      Call DSCAL_(4*lroots**2,-One,SS,1)
      SS(irec(2*lroots,2*lroots))= SS(irec(2*lroots,2*lroots))-One

      End Subroutine mkcipre
