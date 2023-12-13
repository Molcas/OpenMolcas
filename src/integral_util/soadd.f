!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine SOAdd(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,
     &                  iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!***********************************************************************
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      use Constants
      Implicit None
      Integer iBas, jBas, nSOInt, nPrp, lOper,
     &                  iCmp,jCmp,iShell,jShell,iAO,jAO
      Real*8 SOInt(iBas*jBas,nSOInt), PrpInt(nPrp)
      Logical AeqB

      Integer, external:: iPntSO
      Integer i, j, iTri, lSO, j1, i1, j2, j12, i2, iSO1, iSO2, iPnt,
     &        iSO, jSO, Indij, nRow, IndAO1, IndAO2, ip
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!                                                                      *
!***********************************************************************
!                                                                      *

!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In SOAdd:SOInt',' ',SOInt,iBas*jBas,nSOInt)
#endif
!
      lSO = 0
      Do 100 j1 = 0, nIrrep-1
       Do 200 i1 = 1, iCmp
        If (iAOtSO(iAO+i1,j1)<0) Cycle
!
!       Scatter the SO's onto lower rectangular blocks and triangular
!       diagonal blocks.
!
        Do 300 j2 = 0, j1
         j12 = iEor(j1,j2)
         If (iAnd(lOper,2**j12).eq.0) Cycle

         Do 400 i2 = 1, jCmp
          If (iAOtSO(jAO+i2,j2)<0) Cycle
          If (iShell.eq.jShell .and. j1.eq.j2 .and.
     &        i1<i2) Cycle

          lSO = lSO + 1
          iSO1=iAOtSO(iAO+i1,j1)
          iSO2=iAOtSO(jAO+i2,j2)
!
          iPnt = iPntSO(j1,j2,lOper,nbas)
          Do 500 indAO1 = 1, iBas
           Do 600 indAO2 = 1, jBas
            ip = (indAO2-1)*iBas + indAO1
!
!           Move one electron integral.
!
            iSO=iSO1+IndAO1-1
            jSO=iSO2+IndAO2-1

!           Diagonal block. Store only unique elements
            If (j1.eq.j2 .and. iSO1.eq.iSO2 .and.
     &          iSO<jSO) Cycle

            If (j1.eq.j2) Then
!------------Diagonal symmetry block
             Indij=iPnt + iTri(iSO,jSO)
            Else
!------------Off-diagonal symmetry block j1>j2
             nRow = nBas(j1)
             Indij=iPnt + nRow*(jSO-1)*nRow + iSO
            End If

            PrpInt(Indij) = PrpInt(Indij) + SOInt(ip,lSO)
!
 600       Continue
 500      Continue
!
 400     Continue
 300    Continue
!
 200   Continue
 100  Continue
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
        Call Unused_logical(AeqB)
      End If
      End
