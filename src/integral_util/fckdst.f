************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      Subroutine FckDst(TwoHam,nDens,Fij,iBas,jBas,iCmp,jCmp,
     &                  ikop1,ikop2,Irrep,
     &                  Shij,iAO1,iAO2,iAOst1,iAOst2,fact)
      use Symmetry_Info, only: iChTbl
      use SOAO_Info, only: iAOtSO, nSOInf
      Implicit Real*8 (a-h,o-z)
      integer jirr(0:7)
*
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
*
      Real*8 Fij(0:iBas-1,0:jBas-1,iCmp,jCmp),TwoHam(nDens)
      Integer iPnt(0:7)
      Logical Shij
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      iChO=iOper(Irrep)
      If (iChO.eq.0) Then
*
         iPntij = 0
         Do iIrrep = 0, nIrrep-1
            iPnt(iIrrep) = iPntij
            ipntij = ipntij+nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
*
*-----Distribute contributions from the intermediate skeleton
*     Fock matrix onto the symmetry adapted Fock matrix.
*
         iiR = NrOpr(iEor(ikOp1,ikOp2))
         Do i1 = 1, iCmp
          Do i2 = 1, jCmp
           Do iIrrep = 0, nIrrep-1
            If (Max(iAO1+i1,iAO2+i2)>nSOInf) Then
               Write (6,*) 'Fckdst: Max(iSO.jSO)>nSOInf (1)'
               Call Abend()
            End If
            iSO=iAOtSO(iAO1+i1,iIrrep)+iAOst1
            jSO=iAOtSO(iAO2+i2,iIrrep)+iAOst2
            XR = DBLE(iChTbl(iIrrep,iiR))
*
            If (iSO<0 .or. jSO<0) Cycle
*
            ipntij = iPnt(iIrrep)
            Do jAOj = 0, jBas-1
              Do iAOi = 0, iBas-1
                  Fac = XR
                  If (Shij .and. i1.eq.i2 .and.
     &                iAOi+iAOst1.eq.jAOj+iAOst2) Fac = Two*XR
                  jSOj = jSO + jAOj
                  iSOi = iSO + iAOi
                  ipFij = ipntij + iTri(iSOi,jSOj)
                  TwoHam(ipFij) = TwoHam(ipFij)
     &                          + Fact*Fac*Fij(iAOi,jAOj,i1,i2)
              End Do
            End Do

          End Do   ! iIrrep
         End Do    ! i2
        End Do     ! i1
*
      Else
*
         Do iIrrep=0,nIrrep-1
            jIrr(iIrrep)=
     &      NrOpr(iEOr(iOper(iIrrep),iChO))
         End Do
*
         iPntij = 0
         iPnt(0:nIrrep-1)=-1
         Do iIrrep = 0, nIrrep-1
            jIrrep=jIrr(iIrrep)
            If (iIrrep.gt.jIrrep) Then
               iPnt(iIrrep) = iPntij
               iPntij = iPntij+nBas(jIrrep)*nBas(iIrrep)
             End If
         End Do
*
*-----Distribute contributions from the intermediate skeleton
*     Fock matrix onto the symmetry adapted Fock matrix.
*
         l1=NrOpr(ikop1)
         l2=NrOpr(ikop2)
         Do 100 i1 = 1, iCmp
           Do 200 i2 = 1, jCmp
            Do iIrrep = 0, nIrrep-1
               jIrrep=jIrr(iIrrep)
               If (iIrrep.lt.jIrrep) Cycle
               X1 = DBLE(iChTbl(iIrrep,l1))
               X2 = DBLE(iChTbl(jIrrep,l2))
               X3 = DBLE(iChTbl(jIrrep,l1))
               X4 = DBLE(iChTbl(iIrrep,l2))
               If (Max(iAO1+i1,iAO2+i2)>nSOInf) Then
                  Write (6,*) 'Fckdst: Max(iSO.jSO)>nSOInf (2)'
                  Call Abend()
               End If
               iSOi=iAOtSO(iAO1+i1,iIrrep)+iAOst1
               jSOj=iAOtSO(iAO2+i2,jIrrep)+iAOst2
               jSOi=iAOtSO(iAO2+i2,iIrrep)+iAOst2
               iSOj=iAOtSO(iAO1+i1,jIrrep)+iAOst1
               If (iSOi>-1.and.jSOj>-1) Then
                  iPntij = iPnt(iIrrep)
                  Do jAO = 0, jBas-1
                     Do iAO = 0, iBas-1
                        Fac = X1*X2
                        ipF = iPntij
     &                      + (jSOj+jAO-1)*nBas(iIrrep)
     &                      + iSOi + iAO
                        TwoHam(ipF) = TwoHam(ipF)
     &                              + Fact*Fac*Fij(iAO,jAO,i1,i2)
                     End Do
                  End Do
               End If
               If (jSOi>-1.and.iSOj>-1) Then
                  ipntij = iPnt(iIrrep)
                  Do jAO = 0, jBas-1
                     Do iAO = 0, iBas-1
                        Fac = X3*X4
                        ipF = ipntij
     &                      + nBas(iIrrep)*(iSOj+iAO-1)
     &                      + jSOi + jAO
                        TwoHam(ipF) = TwoHam(ipF)
     &                              + Fact*Fac*Fij(iAO,jAO,i1,i2)
                     End Do
                  End Do
               End If
             End Do
 200       Continue
 100     Continue
*
      End If
*
      Return
      End
