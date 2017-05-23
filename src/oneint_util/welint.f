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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine WelInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: to compute the Pauli repulsion integrals with the            *
*         Gauss-Hermite quadrature.                                    *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Rowel                                                   *
*              Traxyz                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. October '92.                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "wldata.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Integer iStabM(0:nStabM-1), lOper(nComp), iChO(nComp)
*
      iRout = 122
      iPrint = nPrint(iRout)
*     iQ = 1
      Call qEnter('WelInt')
      If (iPrint.ge.59) Then
         Write (6,*) ' In WelInt'
         Write (6,*) ' r0, ExpB=',r0,ExpB
         Write (6,*) ' la,lb=',la,lb
      End If
*
      k = la + lb
      jsum = 1
      Do 10 i = 1, k
         jsum = jsum + 3**i
 10   Continue
*
      ip = 1
      ipGri = ip
      ip = ip + nZeta*jsum
      ipGrin= ip
      ip = ip + nZeta*(k+1)*(k/2+1)*(k/4+1)
      iPxyz = ip
      ip = ip + nZeta
      If (ip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2, 'WelInt:  ip-1.gt.nZeta*nArr(pos.1)')
         Write (6,*) ip-1,'>',nZeta*nArr
         Call Abend()
      End If
*
      Call Rowel(nZeta,r0,expB,k,Zeta,P,Array(iPxyz),Array(ipGri),
     &           Array(ipGrin),jsum)
      ip = ip - nZeta
      ip = ip - nZeta*(k+1)*(k/2+1)*(k/4+1)
*
      ipA = ip
      ip = ip + nZeta*9
      ipScr = ip
      ip = ip + nZeta*3**k
      If (ip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2, 'WelInt:  ip-1.gt.nZeta*nArr(pos.2)')
         Write (6,*) ip-1,'>',nZeta*nArr
         Call Abend()
      End If
*
*-----Transform each block to the global coordinate system
*
      iOff = ipgri + nZeta
      Do 100 ik = 1, k
         If (ik.eq.1) Call SetUpA(nZeta,Array(ipA),P)
         Call Traxyz(nZeta,ik,Array(iOff),Array(ipScr),Array(ipA))
         iOff = iOff + nZeta*3**ik
 100  Continue
      If (iPrint.ge.99) Call RecPrt(' In WelInt: Array(ipGri)',' ',
     &   Array(ipGri),nZeta,jSum)
      ip = ip - nZeta*3**k
      ip = ip - nZeta*9
*
      ip1 = ip
      ip = ip + nZeta
      ip2 = ip
      ip = ip + nZeta
      ip3 = ip
      ip = ip + nZeta
      ip4 = ip
      ip = ip + nZeta
      ip5 = ip
      ip = ip + nZeta
      If (ip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2, 'WelInt:  ip-1.gt.nZeta*nArr(pos.3)')
         Write (6,*) ip-1,'>',nZeta*nArr
         Call Abend()
      End If
      Call TraPAB(nZeta,la,lb,Final,Array(ipgri),jSum,rKappa,Array(ip1),
     &            Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,RB,P)
      ip = ip - nZeta*5
      ip = ip - nZeta*jsum
*
      Call qExit('WelInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_integer(nHer)
         Call Unused_real_array(Ccoor)
         Call Unused_real_array(nOrdOp)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_integer_array(iStabM)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
