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
* Copyright (C) 1990,1992,1994,1996, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Tcrtnc(Coef1,n1,m1,Coef2,n2,m2,
     &                  Coef3,n3,m3,Coef4,n4,m4,
     &                  ACInt,mabcd,Scrtch,nScr,ACOut,
     &                  IndZet,lZeta,IndEta,lEta)
************************************************************************
*                                                                      *
* Object: to transform the integrals from primitives to contracted     *
*         basis functions. The subroutine will do both complete and    *
*         incomplete transformations.                                  *
*                                                                      *
*         Observe that ACInt and ACOut may overlap!!!!                 *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to back transformation, January '92.            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
*
*---- Cache size
*
#include "lCache.fh"
      Real*8 Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4),
     &       ACInt(m1*m2*m3*m4,mabcd), Scrtch(nScr),
     &       ACOut(lZeta*lEta,mabcd)
      Integer IndZet(lZeta), IndEta(lEta)
*
      iRout = 18
      iPrint = nPrint(iRout)
      iPrint = 000000 !yma

*
!#ifdef _DEBUGPRINT_  !yma
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(AB|CD)',ACInt,
     &                                m1*m2*m3*m4*mabcd)
      If (iPrint.ge.99) Then
         Call RecPrt(' In Tcrtnc: P(ab|cd)',' ',ACInt,mabcd,m1*m2*m3*m4)
         Call RecPrt(' Coef1',' ',Coef1,n1,m1)
         Call RecPrt(' Coef2',' ',Coef2,n2,m2)
         Call RecPrt(' Coef3',' ',Coef3,n3,m3)
         Call RecPrt(' Coef4',' ',Coef4,n4,m4)
         Write (6,*) n1, n2, n3, n4
      End If
!#endif
*
*---- Reduce for contraction matrix
      nCache = (3*lCache)/4 - n3*m3 - n4*m4
      lsize= m3*m4 + m3*n4
      nVec = m1*m2*mabcd
      IncVec = Min(Max(1,nCache/lsize),nVec)
      ipA3=1
      ipA2=ipA3 + nVec*lEta
*
*define _CHECK_
#ifdef _CHECK_
      If (nVec*lEta+n4*m3*IncVec.gt.nScr) Then
         Write (6,*) 'Tcrtnc: Memory failure 1'
         Write (6,*) 'n4*IncVec*m3(A2)=',n4*IncVec*m3
         Write (6,*) 'nVec*lEta(A3)=',nVec*lEta
         Write (6,*) 'n4,IndVec,m3=',n4,IndVec,m3
         Write (6,*) 'nVec,lEta=',nVec,lEta
         Write (6,*) 'nScr=',nScr
         Call Abend()
      End If
#endif

      Call TncHlf(Coef3,m3,n3,Coef4,m4,n4,lEta,nVec,
     &            IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndEta)
*
      nCache = (3*lCache)/4 - n1*m1 - n2*m2
      lsize = m1*m2 + m1*n2
      nVec = mabcd*lEta
      IncVec = Min(Max(1,nCache/lsize),nVec)
*
#ifdef _CHECK_
      If (nVec*m1*m2+n2*IncVec*m1.gt.nScr) Then
         Write (6,*) 'Tcrtnc: Memory failure 2'
         Write (6,*) 'nVec*m1*m2(A1)=',nVec*m1*m2
         Write (6,*) 'n2*IncVec*m1(A2)=',n2*IncVec*m1
         Write (6,*) 'nVec,m1,m2=',nVec,m1,m2
         Write (6,*) 'n2,IncVec,m1=',n2,IncVec,m1
         Write (6,*) 'nScr=',nScr
         Call Abend()
      End If
#endif
*
      Call TncHlf(Coef1,m1,n1,Coef2,m2,n2,lZeta,nVec,
     &            IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndZet)
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.59)
     &  Call RecPrt(' In Tcrtnc: P(ab|cd) ',' ',ACOut,mabcd,lZeta*lEta)
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(ab|cd)',ACOut,
     &                                lZeta*lEta*mabcd)
#endif
*
*     Call GetMem('Tcrtnc','CHECK','REAL',iDum,iDum)
      Return
      End
      Subroutine Tnchlf(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,
     &                  lZeta,nVec,IncVec,A1,A2,A3,Indij)
************************************************************************
*                                                                      *
* Object: to do a half transformation. The loop over the two matrix-   *
*         matrix multiplications is segmented such that the end of the *
*         intermediate matrix will not push the start of the same out  *
*         from the cache.                                              *
*                                                                      *
* Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*                                                                      *
*             Modified to decontraction May 1996, by R. Lindh          *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2),
     &       A1(nVec,nCntr1,nCntr2), A2(nPrm2,IncVec,nCntr1),
     &       A3(lZeta,nVec)
      Logical Seg1, Seg2
      Integer Indij(lZeta)
*
*-----Check if the basis set is segmented
*
      Seg1=.False.
      Do iPrm1 = nPrm1, 1, -1
         Do iCntr1 = nCntr1, 1, -1
            If (Coeff1(iPrm1,iCntr1).eq.Zero) Then
               Seg1=.True.
               Go To 10
            End If
         End Do
      End Do
 10   Continue
*
      Seg2=.False.
      Do iPrm2 = nPrm2, 1, -1
         Do iCntr2 = nCntr2, 1, -1
            If (Coeff2(iPrm2,iCntr2).eq.Zero) Then
               Seg2=.True.
               Go To 20
            End If
         End Do
      End Do
 20   Continue
*
*-----Set output matrix to zero
*
      Call FZero(A3,nVec*lZeta)
*
*-----Loop sectioning
*
      Do iiVec = 1, nVec, IncVec
         mVec = Min(IncVec,nVec-iiVec+1)
*--------Set intermediate matrix to zero
         Call FZero(A2,nPrm2*IncVec*nCntr1)
*
         If (Seg2) Then
*
*-----First quarter transformation, (x,AB) -> (b,x,a)
*
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
*-----------Check for zero due to segmented basis
            If (Abs(Coeff2(iPrm2,iCntr2)).gt.Zero) Then
               Do iCntr1 = 1, nCntr1
                  Do iVec = 1, mVec
                     A2(iPrm2,iVec,iCntr1) = A2(iPrm2,iVec,iCntr1) +
     &                 Coeff2(iPrm2,iCntr2)
     &                 *A1(iVec+iiVec-1,iCntr1,iCntr2)
                  End Do
               End Do
            End If
         End Do
      End Do
*
         Else    ! Seg2
*
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
            Do iCntr1 = 1, nCntr1
               Do iVec = 1, mVec
                  A2(iPrm2,iVec,iCntr1) = A2(iPrm2,iVec,iCntr1) +
     &              Coeff2(iPrm2,iCntr2)
     &              *A1(iVec+iiVec-1,iCntr1,iCntr2)
               End Do
            End Do
         End Do
      End Do
*
         End If   ! Seg2
*
         If (Seg1) Then
*
*-----Second quarter transformation
*
         Do iCntr1 = 1, nCntr1
            Do iZeta = 1, lZeta
               iPrm2=(Indij(iZeta)-1)/nPrm1 + 1
               iPrm1=Indij(iZeta)-(iPrm2-1)*nPrm1
*--------------Check for zero due to segmented basis
               If (Abs(Coeff1(iPrm1,iCntr1)).gt.Zero) Then
                  ijZeta=(iPrm2-1)*nPrm1 + iPrm1
                  Do iVec = iiVec, iiVec+mVec-1
                     A3(iZeta,iVec) = A3(iZeta,iVec) +
     &                 Coeff1(iPrm1,iCntr1)*
     &                 A2(iPrm2,iVec-iiVec+1,iCntr1)
                  End Do ! iVec
               End If
            End Do       ! iZeta
         End Do          ! iCntr1
*
         Else
*
*-----Second quarter transformation
*
         Do iCntr1 = 1, nCntr1
            Do iZeta = 1, lZeta
               iPrm2=(Indij(iZeta)-1)/nPrm1 + 1
               iPrm1=Indij(iZeta)-(iPrm2-1)*nPrm1
               Do iVec = iiVec, iiVec+mVec-1
                  A3(iZeta,iVec) = A3(iZeta,iVec) +
     &              Coeff1(iPrm1,iCntr1)*
     &              A2(iPrm2,iVec-iiVec+1,iCntr1)
               End Do ! iVec
            End Do    ! iZeta
         End Do       ! iCntr1
*
      End If
*
*-----End of loop sectioning
*
      End Do
*
      Return
      End
