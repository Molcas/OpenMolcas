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
* Copyright (C) 1992,1995, Roland Lindh                                *
************************************************************************
      SubRoutine WelGrd(
#define _CALLING_
#include "grd_interface.fh"
     &                 )
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
*                                                                      *
*             Modified to gradients, April '95. R. Lindh               *
************************************************************************
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "wldata.fh"
#include "print.fh"
#include "disp.fh"

#include "grd_interface.fh"

*
*     Statement function for Cartesian index
*
      nElem(i)=(i+1)*(i+2)/2
*
      iRout = 122
      iPrint = nPrint(iRout)
*     iQ = 1
      If (iPrint.ge.59) Then
         Write (6,*) ' In WelGrd'
         Write (6,*) ' r0, ExpB=',r0,ExpB
         Write (6,*) ' la,lb=',la,lb
         Write (6,*) '  A=',A
         Write (6,*) ' RB=',RB
      End If
*
      k = la + lb + 1
      jsump= 1
      Do i = 1, k
         jsump= jsump+ 3**i
      End Do
      k0 = Max(la + lb - 1,0)
      jsumm= 1
      Do i = 1, k0
         jsumm= jsumm+ 3**i
      End Do
*
      ip = 1
      ipGri = ip
      ip = ip + nZeta*jsump
      ipTGri = ip
      ip = ip + nZeta*jsump
      ipGrin= ip
      ip = ip + nZeta*(k+1)*(k/2+1)*(k/4+1)
      iPxyz = ip
      ip = ip + nZeta
      If (ip-1.gt.nZeta*nArr) Then
         Write (6,*) ' ip-1.gt.nZeta*nArr(pos.1)'
         Write (6,*) ip-1,'>',nZeta*nArr
         Call ErrTra
         Call Abend()
      End If
*
      Call Rowel(nZeta,r0,expB,k,Zeta,P,Array(iPxyz),Array(ipGri),
     &           Array(ipGrin),jsump)
      ip = ip - nZeta
      ip = ip - nZeta*(k+1)*(k/2+1)*(k/4+1)
      If (iPrint.ge.99) Call RecPrt(' In WelInt: Array(ipGri)l',' ',
     &   Array(ipGri),nZeta,jSumP)
*
      ipAMx = ip
      ip = ip + nZeta*9
      ipScr = ip
      ip = ip + nZeta*3**k
      If (ip-1.gt.nZeta*nArr) Then
         Write (6,*) ' ip-1.gt.nZeta*nArr(pos.2)'
         Write (6,*) ip-1,'>',nZeta*nArr
         Call ErrTra
         Call Abend()
      End If
*
*-----Transform each block to the global coordinate system
*
      iOff = ipgri + nZeta
      Do 100 ik = 1, k
         If (ik.eq.1) Call SetUpA(nZeta,Array(ipAMx),P)
         Call Traxyz(nZeta,ik,Array(iOff),Array(ipScr),Array(ipAMx))
         iOff = iOff + nZeta*3**ik
 100  Continue
      If (iPrint.ge.99) Call RecPrt(' In WelInt: Array(ipGri)g',' ',
     &   Array(ipGri),nZeta,jSumP)
      call dcopy_(nZeta*jsump,Array(ipGri),1,Array(ipTGri),1)
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
         Write (6,*) ' ip-1.gt.nZeta*nArr(pos.3)'
         Write (6,*) ip-1,'>',nZeta*nArr
         Call ErrTra
         Call Abend()
      End If
*
*----- Compute <a|O|b+1>
*
      ip0p=ip
      ip  =ip + nZeta*nElem(la)*nElem(lb+1)
      Call TraPAB(nZeta,la,lb+1,Array(ip0p),Array(ipgri),jSump,rKappa,
     &     Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,
     &     RB,P)
*
*----- Compute <a|O|b-1>
*
      If (lb.ge.1) Then
         ip0m=ip
         ip  =ip + nZeta*nElem(la)*nElem(lb-1)
         call dcopy_(nZeta*jsumm,Array(ipTGri),1,Array(ipGri),1)
         Call TraPAB(nZeta,la,lb-1,Array(ip0m),Array(ipgri),jSumm,
     &               rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),
     &               Array(ip5),A,RB,P)
      Else
         ip0m=1
      End If
*
*----- Compute <a+1|O|b>
*
      ipp0=ip
      ip  =ip + nZeta*nElem(la+1)*nElem(lb)
      call dcopy_(nZeta*jsump,Array(ipTGri),1,Array(ipGri),1)
      Call TraPAB(nZeta,la+1,lb,Array(ipp0),Array(ipgri),jSump,rKappa,
     &     Array(ip1),Array(ip2),Array(ip3),Array(ip4),Array(ip5),A,
     &     RB,P)
*
*----- Compute <a-1|O|b>
*
      If (la.ge.1) Then
         ipm0=ip
         ip  =ip + nZeta*nElem(la-1)*nElem(lb)
         call dcopy_(nZeta*jsumm,Array(ipTGri),1,Array(ipGri),1)
         Call TraPAB(nZeta,la-1,lb,Array(ipm0),Array(ipgri),jSumm,
     &               rKappa,Array(ip1),Array(ip2),Array(ip3),Array(ip4),
     &               Array(ip5),A,RB,P)
      Else
         ipm0=1
      End If
*
      ipAlph = ip
      ip = ip + nZeta
      ipBeta = ip
      ip = ip + nZeta
*
      jp = ipAlph
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(jp),1)
         jp = jp + nAlpha
      End Do
      jp = ipBeta
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(jp),nAlpha)
         jp = jp + 1
      End Do
*
*---- Assemble the derivative integrals and distribute contributions.
*
      Call CmbnW1(Array(ipp0),Array(ipm0),Array(ip0p),Array(ip0m),
     &            nZeta,la,lb,Zeta,rKappa,Final,Array(ipAlph),
     &            Array(ipBeta),Grad,nGrad,DAO,
     &            IfGrad,IndGrd,dc(mdc)%nStab,dc(ndc)%nStab,kOp)

*
      ip = ip - 5*nZeta
      ip = ip - 2*nZeta
      If (la.ge.1) ip = ip - nZeta*nElem(la-1)*nElem(lb)
      ip = ip - nZeta*nElem(la)*nElem(lb+1)
      If (lb.ge.1) ip = ip - nZeta*nElem(la)*nElem(lb-1)
      ip = ip - nZeta*nElem(la+1)*nElem(lb)
      ip = ip - 2*nZeta*jsump
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer(nHer)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iStabM)
      End If
      End
