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
************************************************************************
      SubRoutine Cntrct(First,
     &                  Coef1,n1,m1,Coef2,n2,m2,
     &                  Coef3,n3,m3,Coef4,n4,m4,
     &                  ACInt,mabMin,mabMax,mcdMin,mcdMax,
     &                  Scrtch,nScrtch,ACOut,
     &                  IndZet,lZeta,IndEta,lEta)
************************************************************************
*                                                                      *
* Object: to transform the integrals from primitives to contracted     *
*         basis functions. The subroutine will do both complete and    *
*         incomplete transformations.                                  *
*                                                                      *
* Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
*
*-----Cache size
*
#include "lCache.fh"
      Real*8 Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4),
     &       ACInt(n1*n2*n3*n4,(mabMax-mabMin+1)*(mcdMax-mcdMin+1)),
     &       Scrtch(nScrtch),
     &       ACOut((mabMax-mabMin+1)*(mcdMax-mcdMin+1),m1*m2*m3*m4)
      Logical First
      Integer IndZet(lZeta), IndEta(lEta)
*
      iRout = 18
      iPrint = nPrint(iRout)
*
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
         Call RecPrt('Cntrct: Coef1',' ',Coef1,n1,m1)
         Call RecPrt('Cntrct: Coef2',' ',Coef2,n2,m2)
         Call RecPrt('Cntrct: Coef3',' ',Coef3,n3,m3)
         Call RecPrt('Cntrct: Coef4',' ',Coef4,n4,m4)
         Call RecPrt('Cntrct: [a0|c0]',' ',ACInt,lZeta,lEta*mabcd)
         Write (6,*) 'IndZet=',IndZet
         Write (6,*) 'IndEta=',IndEta
      If (.not.First)
     &   Call RecPrt(' In Cntrct: Partial (a0|c0)',' ',
     &               ACOut,mabcd,m1*m2*m3*m4)
#endif
*     The idea here is to make the transformation in subblocks
*     (size=IncVec) to minimize cache faults. We split the range of
*     the compound index such that the contraction coefficients
*     (Coef1{n1,m1} & Coef2{n2,m2}), the first quater transformed block
*     {n2,IncVec}, for a fixed m1 index, and the half transformed block
*     {IncVec,m1,m2} fit into a fixed cache size.
*
*-----Reduce for contraction matrices and 3/4th
      nCache_ = (3*lCache)/4 - n1*m1 - n2*m2
*     Compute the size of the first quarter and half transformed block.
      lsize =  n1*n2 + n2*m1
*     The length of the compound index
      nVec = lEta*mabcd
*     Compute the size of the increment of which we will run the
*     compound index. It is the largest of 1 or the
      IncVec = Min(Max(1,nCache_/lsize),nVec)
*     Pointer to the full block of half transformed integrals
*     {nVec*m1*m2}
      ipA3 = 1
*     Pointer to the first quater transformed block {n1*IncVec}
      ipA2 = ipA3 + nVec*m1*m2
*define _CHECK_
#ifdef _CHECK_
      If (nVec*m1*m2+n2*IncVec.gt.nScrtch) Then
         Write (6,*) 'Cntrct: Memory failure 1'
         Write (6,*) 'nVec*m1*m2(A3)=',nVec*m1*m2
         Write (6,*) 'n2*IncVec(A2)=',n2*IncVec
         Write (6,*) 'n2,IncVec=',n2,IncVec
         Write (6,*) 'nVec,lsize=',nVec,lSize
         Write (6,*) 'n1,n2,m1  =',n1, n2, m1
         Write (6,*) 'nScrtch=',nScrtch
         Call Abend()
      End If
#endif
*
      Call CntHlf(Coef1,m1,n1,Coef2,m2,n2,lZeta,nVec,
     &            .True.,IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),
     &            IndZet)
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Call RecPrt('Halftransformed',' ',
     &                       Scrtch(ipA3),nVec,m1*m2)
#endif
*
      nCache_ = (3*lCache)/4 - n3*m3 - n4*m4
      lsize = n3*n4 + n4*m3
      nVec = mabcd*m1*m2
      IncVec = Min(Max(1,nCache_/lsize),nVec)
#ifdef _CHECK_
      If (nVec*m3*m4+n4*IncVec.gt.nScrtch) Then
         Write (6,*) 'Cntrct: Memory failure 2'
         Call Abend()
      End If
#endif
*
      Call CntHlf(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,
     &            First,IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,
     &            IndEta)
      First = .False.
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In Cntrct: (a0|c0) ',' ',
     &            ACOut,mabcd,m1*m2*m3*m4)
#endif
*
      Return
      End
