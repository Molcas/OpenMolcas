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
      Subroutine Fluctuating(AInv,nAtoms,Lambda,dQ,nij,nPert,
     &                       iANr,rMP,nElem,EC,Alpha)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "constants.fh"
      Real*8 AInv(nAtoms,nAtoms), Lambda(nAtoms), dQ(nAtoms), EC(3,nij),
     &       A(3), B(3), rMP(nij,0:nElem-1,0:nPert-1)
      Integer iANr(nAtoms)
*
      Do iPert = 1, 6
         Do iAtom = 1, nAtoms
            ii=iAtom*(iAtom+1)/2
            dQ(iAtom)=rMP(ii,0,0)-rMP(ii,0,iPert)
         End Do
c         Call RecPrt('dQ',' ',dQ,1,nAtoms)
         Call DGEMM_('N','N',
     &               nAtoms,1,nAtoms,
     &               1.0d0,AInv,nAtoms,
     &               dQ,nAtoms,
     &               0.0d0,Lambda,nAtoms)
c         Call RecPrt('Lambda',' ',Lambda,1,nAtoms)
*
         Do iAtom = 1, nAtoms
            R_BS_i=Bragg_Slater(iANr(iAtom))
            ii = iAtom*(iAtom+1)/2
            call dcopy_(3,EC(1,ii),1,A,1)
*
            Do jAtom = 1, iAtom-1
               R_BS_j=Bragg_Slater(iANr(jAtom))
               jj = jAtom*(jAtom+1)/2
               call dcopy_(3,EC(1,jj),1,B,1)
               rij2=(A(1)-B(1))**2 +(A(2)-B(2))**2 +(A(3)-B(3))**2
               ri=Lambda(iAtom)
               rj=Lambda(jAtom)
               rij02=((R_BS_i+R_BS_j))**2
               ij=iAtom*(iAtom-1)/2+jAtom
C              Write (6,*) ij,ri,rj,rij02
               rMP(ij,0,iPert)=-(ri-rj)*Exp(-Alpha*(rij2/rij02))/Two
C               Write (*,*) 'RMP',iAtom,jAtom,rMP(ij,0,iPert)
            End Do
*
         End Do
*
      End Do
c      Call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
*
      Return
      End
