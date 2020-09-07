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
      Subroutine ProjSym(nAtoms,nCent,Ind,nStab,jStab,A,
     &                   iDCRs,B,Smmtrc,nDim,Print,dB,
     &                   mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,
     &                   nB_Tot,ndB_Tot,Proc_dB,nqB,nB,iq,
     &                   rMult)
      use Phase_Info
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
*
#include "real.fh"
#include "print.fh"
      Real*8 Tx(3,MxAtom), A(3,nCent), B(3,nCent),
     &       dB(3,nCent,3,nCent), BM(nB_Tot), dBM(ndB_Tot)
      Integer   Ind(nCent), nStab(nAtoms), jStab(0:7,nAtoms),
     &          iDCRs(nCent), iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
      Logical Smmtrc(3,nAtoms), Print, Proc_dB
*
      If (Print) Then
         Call RecPrt('B',' ',B,3,nCent)
         Call RecPrt('dB',' ',dB,3*nCent,3*nCent)
         Write (6,*) iDCRs
      End If
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
         Do ixyz = 1, 3
            Tx(ixyz,i)=DBLE(iPhase(ixyz,iDCRs(i)))*Tx(ixyz,i)
         End Do
      End Do
*
*---- Create BqR
*
      nq = 0
      Do i = 1, nCent
         Do ixyz = 1, 3
            If (Smmtrc(ixyz,Ind(i))) Then
               iDim = 0
               Do j = 1, Ind(i)
                  jxyz_Max=3
                  If (j.eq.Ind(i)) jxyz_Max=ixyz
                  Do jxyz = 1, jxyz_Max
                     If (Smmtrc(jxyz,j)) iDim = iDim + 1
                  End Do
               End Do
               mB_Tot = mB_Tot + 1
               nq = nq + 1
               BM(mB_Tot) = Tx(ixyz,i)*B(ixyz,i)
               iBM(mB_Tot) = iDim
            End If
         End Do
      End Do
      nqB(iq) = nq
*
*---- Create dBqR
*
      If (.Not.Proc_dB) Go To 99
      Do i = 1, nCent
         Do ixyz = 1, 3
            If (Smmtrc(ixyz,Ind(i))) Then
               iDim = 0
               Do j = 1, Ind(i)
                  jxyz_Max=3
                  If (j.eq.Ind(i)) jxyz_Max=ixyz
                  Do jxyz = 1, jxyz_Max
                     If (Smmtrc(jxyz,j)) iDim = iDim + 1
                  End Do
               End Do
*
               Do j = 1, nCent
                  Do jxyz = 1, 3
                     If (Smmtrc(jxyz,Ind(j))) Then
*
                        jDim = 0
                        Do k = 1, Ind(j)
                           kxyz_Max=3
                           If (k.eq.Ind(j)) kxyz_Max=jxyz
                           Do kxyz = 1, kxyz_Max
                              If (Smmtrc(kxyz,k)) jDim = jDim + 1
                           End Do
                        End Do
*
                        mdB_Tot = mdB_Tot + 1
                        dBM(mdB_Tot) = rMult*
     &                      Tx(ixyz,i)*dB(ixyz,i,jxyz,j)*Tx(jxyz,j)
                        idBM(1,mdB_Tot) = iDim
                        idBM(2,mdB_Tot) = jDim
*
                     End If
                  End Do
            End Do

            End If
         End Do
      End Do
 99   Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nDim)
      End
      Subroutine ProjSym2(nAtoms,nCent,Ind,nStab,jStab,A,
     &                   iDCRs,B,BqR,Smmtrc,Print,dB,dBqR)
      use Phase_Info
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
*
#include "real.fh"
#include "print.fh"
      Real*8 Tx(3,MxAtom), A(3,nCent), B(3,nCent), BqR(3,nAtoms),
     &       dB(3,nCent,3,nCent), dBqR(3,nAtoms,3,nAtoms)
      Integer   Ind(nCent), nStab(nAtoms), jStab(0:7,nAtoms),
     &          iDCRs(nCent)
      Logical Smmtrc(3,nAtoms), Print
*
      If (Print) Then
         Call RecPrt('B',' ',B,3,nCent)
         Call RecPrt('dB',' ',dB,3*nCent,3*nCent)
         Write (6,*) iDCRs
      End If
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
         Do ixyz = 1, 3
            Tx(ixyz,i)=DBLE(iPhase(ixyz,iDCRs(i)))*Tx(ixyz,i)
         End Do
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
      If (Print) Call RecPrt('BqR',' ',BqR,1,3*nAtoms)
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
      If (Print) Call RecPrt('dBqR',' ',dBqR,3*nAtoms,3*nAtoms)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical_array(Smmtrc)
      End
