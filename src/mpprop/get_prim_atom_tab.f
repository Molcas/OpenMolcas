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
      Subroutine Get_Prim_Atom_Tab(nAtoms,nPrim,
     &Work,CENTX,CENTY,CENTZ)
      Implicit Real*8 (a-h,o-z)

#include "MolProp.fh"
      Dimension Work(3*nAtoms)
! Begin EB
! EB     Dimension CENTX(nPrim),CENTY(nPrim),CENTZ(nPrim)
      Dimension CENTX(nPrim*(nPrim+1)/2),
     &          CENTY(nPrim*(nPrim+1)/2),
     &          CENTZ(nPrim*(nPrim+1)/2)
! End EB
*
*---- Get the primitiv basis that belongs to a specific atom
*
      Do i=1,nAtoms
         nAtomPBas(i) = 0
         Do j=1,nPrim
            If((ABS(Work((i-1)*3+1)-CENTX(j*(j-1)/2+j)).le.1E-10).and.
     &         (ABS(Work((i-1)*3+2)-CENTY(j*(j-1)/2+j)).le.1E-10).and.
     &         (ABS(Work((i-1)*3+3)-CENTZ(j*(j-1)/2+j)).le.1E-10))Then
               nAtomPBas(i) = nAtomPBas(i) + 1
               iAtPrTab(nAtomPBas(i),i) = j
!               Write(*,*)'Atom tab',i,j,CENTX(j*(j-1)/2+j),
!     &         CENTY(j*(j-1)/2+j),CENTZ(j*(j-1)/2+j),nAtomPBas(i)
            EndIf
         EndDo
      EndDo

      Return
      End
