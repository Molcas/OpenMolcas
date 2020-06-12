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
* Copyright (C) 1991,2001, Roland Lindh                                *
************************************************************************
      SubRoutine PCMInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,CCoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              mHrr                                                    *
*              DCR                                                     *
*              Rys                                                     *
*              Hrr                                                     *
*              DaXpY   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, January '91                             *
*                                                                      *
*             Modified to PCM-integrals, by RL June '01, Napoli, Italy.*
************************************************************************
      use PCM_arrays
      Implicit Real*8 (A-H,O-Z)
*     Used for normal nuclear attraction integrals
      External TNAI, Fake, XCff2D, XRys2D
*     Used for finite nuclei
      External TERI, ModU2, vCff2D, vRys2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,nComp),
     &       Array(nZeta*nArr)
      Integer iStabM(0:nStabM-1), lOper(nComp)
*-----Local arrys
      Real*8 C(3), TC(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga(4), iDCRT(0:7), iChO(nComp)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
      Dimension jStab_(0:0)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 151
      iPrint = nPrint(iRout)
*     Call qEnter('PCMInt')
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(2*3,Coora,1,Coori,1)
      mabMin = nabSz(Max(la,lb)-1)+1
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mabMax = nabSz(la+lb)
*
*     Compute FLOP's and size of work array which Hrr will use.
*
      Call mHrr(la,lb,nFLOP,nMem)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
*---- The coordinates of the individual tiles are not stabilized by any
*     operator but the unit operator.
*
      nStab_=1
      jStab_=0
*
*     Loop over tiles.
*
      Do iTile = 1, nTiles
         QTessera = Q_Tessera(iTile)
         C(:)=C_Tessera(:,iTile)
         If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
*--------Find the DCR for M and S
*
         Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &            jStab_ ,nStab_,iDCRT,nDCRT)
         Fact = One / DBLE(LmbdT)
*
         If (iPrint.ge.99) Then
            Write (6,*) ' m      =',nStabM
            Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &            ii = 0, nStabM-1)
            Write (6,*) ' s      =',nStab_
            Write (6,'(9A)') '(S)=',ChOper(jStab_)
            Write (6,*) ' LambdaT=',LmbdT
            Write (6,*) ' t      =',nDCRT
            Write (6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),
     &            ii = 0, nDCRT-1)
         End If

*
         Do lDCRT = 0, nDCRT-1
            TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
            TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
            TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
            call dcopy_(3,TC,1,Coora(1,3),1)
            call dcopy_(3,TC,1,Coora(1,4),1)
*
*           Compute integrals with the Rys quadrature.
*
            nT = nZeta
            NoSpecial=.True.
            Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &               [One],[One],1,P,nZeta,
     &               TC,1,rKappa,[One],Coori,Coora,CoorAC,
     &               mabmin,mabmax,0,0,Array,nArr*nZeta,
     &               TNAI,Fake,XCff2D,XRys2D,NoSpecial)
*
*-----------Use the HRR to compute the required primitive integrals.
*
            Call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)
*
*-----------Accumulate contributions to the symmetry adapted operator
*
            nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
            Call SymAdO(Array(ipIn),nZeta,la,lb,nComp,Final,nIC,
     &                  nOp         ,lOper,iChO,-Fact*QTessera)
            If (iPrint.ge.99) Then
               Write (6,*) Fact*QTessera
               Call RecPrt('PCMInt: Array(ipIn)',' ',Array(ipIn),
     &                 nZeta,nElem(la)*nElem(lb)*nComp)
               Call RecPrt('PCMInt: Final',' ',Final,
     &                 nZeta,nElem(la)*nElem(lb)*nIC)
            End If
*
         End Do
      End Do
*
*     Call GetMem(' Exit PCMInt','LIST','REAL',iDum,iDum)
*     Call qExit('PCMInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_real_array(CCoor)
         Call Unused_integer(nOrdOp)
      End If
      End
