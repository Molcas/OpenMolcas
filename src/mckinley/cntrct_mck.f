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
* Copyright (C) 1994, Roland Lindh                                     *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Cntrct_mck(First,
     &                  Coef1,n1,m1,Coef2,n2,m2,
     &                  Coef3,n3,m3,Coef4,n4,m4,
     &                  g1In,nGr,Array,nArr,
     &                  xpre,G1Out,ngr1,nt,
     &                  IndZet,nZeta,lZeta,IndEta,nEta,lEta)
************************************************************************
*                                                                      *
* Object: to transform the integrals from primitives to contracted     *
*         basis functions. The subroutine will do both complete and    *
*         incomplete transformations.                                  *
*                                                                      *
* Called from: RysG2                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CntHlf                                                  *
*              QExit                                                   *
*                                                                      *
* Author:      Roland Lindh, Dept. of Theoretical Chemistry, University*
*              of Lund, SWEDEN.                                        *
*                                                                      *
* Modified by: Anders Bernhardsson for direct implementation of the    *
*              calculation of first order derivatives needed for       *
*              response calculation.                                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
c#include "print.fh"
#include "lCache.fh"
      Real*8 Coef1(n1,m1),Coef2(n2,m2),Coef3(n3,m3),Coef4(n4,m4),
     &       g1In(nT,nGr),
     &       Array(nArr),
     &       g1Out(nGr1),xpre(nt)
      Logical First
      Integer IndZet(nZeta), IndEta(nEta)
*
c     iRout = 18
c     iPrint = nPrint(iRout)
c     Call QEnter('Cntrct')
*
c     If (iPrint.ge.99)
c    &   Call RecPrt(' In Cntrct: ',' ',G1In,nt,nGr)
c     If (iPrint.ge.59 .and. .not.First)
c    &   Call RecPrt(' In Cntrct: Partial (a0|c0)',' ',
c    &               G1Out,nGr,m1*m2*m3*m4)

*-----Cache size is 32 k word (real*8)
*
      Do iabcdg=1,ngr
         Do it=1,nt
           G1In(it,iabcdg)=G1In(it,iabcdg)*xpre(it)
         End Do
      End Do
*
*-----Reduce for contraction matrix
      nCache = (3*lCache)/4 - n1*m1 - n2*m2
      lsize =  n1*n2 + n2*m1
      nVec = lEta*nGr
      IncVec = Min(Max(1,nCache/lsize),nVec)
      ipA3 = 1
      ipA2 = ipA3 + nVec*m1*m2
      ip=ipA2 + n2*IncVec*m1
      If (ip.gt.nArr)       Call Abend
*
      Call CntHlf_mck(Coef1,m1,n1,Coef2,m2,n2,nZeta,lZeta,nVec,
     &                .True.,IncVec,G1In,Array(ipA2),Array(ipA3),IndZet)
*
      nCache = (3*lCache)/4 - n3*m3 - n4*m4
      lsize = n3*n4 + n4*m3
      nVec = nGr*m1*m2
      IncVec = Min(Max(1,nCache/lsize),nVec)
      ip=ipA2 + n4*IncVec*m3
      If (ip.gt.nArr)       Call Abend
*
      Call CntHlf_mck(Coef3,m3,n3,Coef4,m4,n4,nEta,lEta,nVec,
     &                First,IncVec,Array(ipA3),Array(ipA2),G1Out,IndEta)
      First = .False.
*
c     If (iPrint.ge.59)
c    &  Call RecPrt(' In Cntrct:  ',' ',
c    &              ACOut,labcdG,m1*m2*m3*m4)
*
*     Call GetMem('Cntrct','Check','Real',iDum,iDum)
c     Call QExit('Cntrct')
      Return
      End
