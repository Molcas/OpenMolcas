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
* Copyright (C) 1991,1995, Roland Lindh                                *
************************************************************************
      SubRoutine EFInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                 iStabM,nStabM,
     &                 PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of electric field         *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy  (ESSL)                                           *
*              SOS                                                     *
*              DCR                                                     *
*              XRys                                                    *
*              Util1                                                   *
*              DaXpY  (ESSL)                                           *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, January '91                             *
*                                                                      *
* Modified for explicit code, R. Lindh, February '95.                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake, XCff2D, XRys2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7),
     &        iStabO(0:7), lOper(nComp), iChO(nComp)
*---- Local arrays
      Real*8 Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga(4)
      Character*80 Label
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 200
      iPrint = nPrint(iRout)
      Call qEnter('EFInt')
*
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = nOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      mabMin=nabSz(Max(la,lb)-1)+1
      mabMax=nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mcdMin=nabSz(nOrdOp-1)+1
      mcdMax=nabSz(nOrdop)
      lab=(mabMax-mabMin+1)
      kab=nElem(la)*nElem(lb)
      lcd=(mcdMax-mcdMin+1)
      labcd=lab*lcd
*
*---- Compute Flop's and size of work array which HRR will Use.
*
      Call mHRR(la,lb,nFLOP,nMem)
*
*---- Distribute the work array
*
      ip2 = 1
      ip1 = ip2 + nZeta*Max(labcd,lcd*nMem)
      mArr = nArr - Max(labcd,lcd*nMem)
*
*---- Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
         call dcopy_(3, A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
*
*
      Do 102 lDCRT = 0, nDCRT-1
         TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1)
         TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2)
         TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3)
         call dcopy_(3,TC,1,CoorAC(1,2),1)
         call dcopy_(3,TC,1, Coori(1,3),1)
         call dcopy_(3,TC,1, Coori(1,4),1)
*
*------- Compute integrals with the Rys-Gauss quadrature.
*
         nT=nZeta
         NoSpecial=.True.
         Call Rys(iAnga,nT,
     &            Zeta,ZInv,nZeta,[One],[One],1,
     &            P,nZeta,TC,1,
     &            rKappa,[One],Coori,Coori,CoorAC,
     &            mabMin,mabMax,mcdMin,mcdMax,Array(ip1),mArr*nZeta,
     &            TNAI,Fake,XCff2D,XRys2D,NoSpecial)
*
*------- The integrals are now ordered as ijkl,e,f
*
*        a) Change the order to f,ijkl,e
*        b) Unfold e to ab, f,ijkl,ab
*        c) Change the order back to ijkl,ab,f
*
*a)-----
*
         Call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)
*
*b)----- Use the HRR to unfold e to ab
*
         Call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
         ip3=ip2-1+ipIn
*
*c)-----
*
         Call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,Array(ip1),nZeta*kab)
*
*------- Modify to traceless form, the sixth element contains r*r and
*
         If (nOrdOp.eq.2) Then
*        If (.False.) Then
            nzab=nZeta*kab
            iOffxx=ip1
            iOffyy=ip1+nzab*3
            iOffzz=ip1+nzab*5
            ThreeI = One / Three
            Do i = 0, nzab-1
               RR = Array(iOffxx+i)
     &            + Array(iOffyy+i)
     &            + Array(iOffzz+i)
               XX = Two * Array(iOffxx+i) -
     &            Array(iOffyy+i) - Array(iOffzz+i)
               YY = Two * Array(iOffyy+i) -
     &            Array(iOffxx+i) - Array(iOffzz+i)
               ZZ = Two * Array(iOffzz+i) -
     &            Array(iOffxx+i) - Array(iOffyy+i)
               Array(iOffxx+i) = XX * ThreeI
               Array(iOffyy+i) = YY * ThreeI
               Array(iOffzz+i) = RR
            End Do
         End If
*
*        Stored as nZeta,iElem,jElem,iComp
*
         If (iPrint.ge.49) Then
            Write (6,*) ' In EFInt la,lb=',la,lb
            nzab=nZeta*kab
            Do iElem = 1, nElem(la)
               Do jElem = 1, nElem(lb)
                  ij = (jElem-1)*nElem(la) + iElem
                  ip = ip1 + nZeta*(ij-1)
                  Do iComp = 1, nComp
                     Write (Label,'(A,I2,A,I2,A,I2,A)')
     &                     ' Final (',iElem,',',jElem,',',iComp,') '
                     Call RecPrt(Label,' ',Array(ip),nZeta,1)
                     ip = ip + nzab
                  End Do
               End Do
            End Do
         End If
*
*------- Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ip1),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 102  Continue
*     Call GetMem(' Exit EFInt','LIST','REAL',iDum,iDum)
      Call qExit('EFInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
