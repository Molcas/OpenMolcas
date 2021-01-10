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
      Subroutine ProjSym2(nAtoms,nCent,Ind,A,iDCRs,B,BqR,dB,dBqR)
      use Slapaf_Info, only: jStab,nStab
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
*
#include "real.fh"
      Real*8 Tx(3,MxAtom), A(3,nCent), B(3,nCent), BqR(3,nAtoms),
     &       dB(3,nCent,3,nCent), dBqR(3,nAtoms,3,nAtoms), ATemp(3)
      Integer   Ind(nCent), iDCRs(nCent)
*
#ifdef _DEBUGPRINT_
      Call RecPrt('B',' ',B,3,nCent)
      Call RecPrt('dB',' ',dB,3*nCent,3*nCent)
      Write (6,*) iDCRs
#endif
*
*---- Set up the T-matrix
*
*---- Project away nonsymmetric displacements
*
      call dcopy_(3*nCent,[One],0,Tx,1)
      Do i = 1, nCent
         Call NonSym(nStab(Ind(i)),jStab(0,Ind(i)),A(1,i),Tx(1,i))
*
*------- Rotate vector back to the unique center
*
         Call OA(iDCRs(i),Tx(1:3,i),ATemp)
         Tx(:,i)=ATemp(:)
      End Do
*
*---- The T-matrix is now computed. Now create BqR and dBqR.
*
*---- Create BqR
*
      Call FZero(BqR,3*nAtoms)
      Do i = 1, nCent
         Do ixyz = 1, 3
            BqR(ixyz,Ind(i)) = BqR(ixyz,Ind(i))
     &                       + Tx(ixyz,i)*B(ixyz,i)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('BqR',' ',BqR,1,3*nAtoms)
#endif
*
*---- Create dBqR
*
      Call FZero(dBqR,(3*nAtoms)**2)
      Do i = 1, nCent
         Do ixyz = 1, 3
*
            Do j = 1, nCent
               Do jxyz = 1, 3
*
                  dBqR(ixyz,Ind(i),jxyz,Ind(j)) =
     &               dBqR(ixyz,Ind(i),jxyz,Ind(j))
     &             + Tx(ixyz,i)*dB(ixyz,i,jxyz,j)*Tx(jxyz,j)
*
               End Do
            End Do

         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('dBqR',' ',dBqR,3*nAtoms,3*nAtoms)
#endif
*
      Return
      End
