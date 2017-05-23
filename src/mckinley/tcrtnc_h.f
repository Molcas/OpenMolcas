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
      SubRoutine Tcrtnc_h(Coef1,n1,m1,Coef2,n2,m2,
     &                    Coef3,n3,m3,Coef4,n4,m4,
     &                    ACInt,mabcd,Scrtch,nScr,ACOut,
     &                    IndZet,lZeta,IndEta,lEta)
************************************************************************
*                                                                      *
* Object: to transform the integrals from primitives to contracted     *
*         basis functions. The subroutine will do both complete and    *
*         incomplete transformations.                                  *
*                                                                      *
*         Observe that ACInt and ACOut may overlap!!!!                 *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              Trnglr                                                  *
*              DGEMM_  (ESSL)                                          *
*              DGeTMO  (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
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
     &       ACOut(n1*n2*n3*n4,mabcd)
      Integer IndZet(lZeta), IndEta(lEta)
*
      iRout = 18
      iPrint = nPrint(iRout)
*     iPrint=99
*     Call qEnter('Tcrtnc')
*
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(AB|CD)',ACInt,
     &                                m1*m2*m3*m4*mabcd)
      If (iPrint.ge.99) Then
         Call RecPrt(' In Tcrtnc: P(ab|cd)',' ',ACInt,m1*m2,m3*m4*mabcd)
         Call RecPrt(' Coef1',' ',Coef1,n1,m1)
         Call RecPrt(' Coef2',' ',Coef2,n2,m2)
         Call RecPrt(' Coef3',' ',Coef3,n3,m3)
         Call RecPrt(' Coef4',' ',Coef4,n4,m4)
         Write (6,*) n1, n2, n3, n4
      End If
*
*---- Reduce for contraction matrix
      nCache = (3*lCache)/4 - n1*m1 - n2*m2
      lsize= m1*m2 + m2*n1
      nVec = m3*m4*mabcd
      IncVec = Min(Max(1,nCache/lsize),nVec)
      ipA3=1
      nA3=nVec*lZeta         ! This is the same for the second set!
      ipA2=ipA3 + nA3
      nA2=IncVec*n1*m2
*
      Call TncHlf_h(Coef1,m1,n1,Coef2,m2,n2,iDum,lZeta,nVec,
     &              IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndZet)
*
      nCache = (3*lCache)/4 - n3*m3 - n4*m4
      lsize = m3*m4 + m4*n3
      nVec = mabcd*lZeta
      IncVec = Min(Max(1,nCache/lsize),nVec)
*
      lZE=lZeta*lEta
      If (mabcd.ne.1) Then
         Call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,iDum,lEta,nVec,
     &                 IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
         Call DGeTMO(ACOut,mabcd,mabcd,lZE,Scrtch,lZE)
         call dcopy_(mabcd*lZE,Scrtch,1,ACOut,1)
      Else
         Call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,iDum,lEta,nVec,
     &                 IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
      End If
*
      If (iPrint.ge.59)
     &  Call RecPrt(' In Tcrtnc: P(ab|cd) ',' ',ACOut,mabcd,lZE)
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(ab|cd)',ACOut,
     &                                lZE*mabcd)
*
*     Call GetMem('Tcrtnc','CHECK','REAL',iDum,iDum)
*     Call qExit('Tcrtnc')
      Return
      End
      Subroutine Tnchlf_h(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,
     &                    nPrm2,mZeta,lZeta,nVec,IncVec,A1,A2,A3,Indij)
************************************************************************
*                                                                      *
* Object: to do a half transformation. The loop over the two matrix-   *
*         matrix multiplications is segmented such that the end of the *
*         intermediate matrix will not push the start of the same out  *
*         from the cache.                                              *
*                                                                      *
* Called from: Cntrct                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
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
     &       A1(nCntr1,nCntr2,nVec), A2(nCntr2,IncVec*nPrm1),
     &       A3(nVec,lZeta)
      Integer Indij(lZeta)
      Logical Seg1, Seg2
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
         call dcopy_(nCntr2*mVec*nPrm1,Zero,0,A2,1)
*
         If (Seg1) Then
*
*-----First quarter transformation
*
*
      Do iPrm1 = 1, nPrm1
         Do iCntr1 = 1, nCntr1
*-----------Check for zero due to segmented basis
            If (Abs(Coeff1(iPrm1,iCntr1)).gt.Zero) Then
               Do iCntr2 = 1, nCntr2
                  Do iVec = iiVec, iiVec+mVec-1
                     ijVec = mVec*(iPrm1-1) + (iVec-iiVec+1)
                     A2(iCntr2,ijVec) = A2(iCntr2,ijVec) +
     &                 Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iVec)
                  End Do
               End Do
            End If
         End Do
      End Do
*
         Else    ! Seg1
*
*-----First quarter transformation
*
*
      Do iPrm1 = 1, nPrm1
         Do iCntr1 = 1, nCntr1
            Do iCntr2 = 1, nCntr2
               Do iVec = iiVec, iiVec+mVec-1
                  ijVec = mVec*(iPrm1-1) + (iVec-iiVec+1)
                  A2(iCntr2,ijVec) = A2(iCntr2,ijVec) +
     &              Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iVec)
               End Do
            End Do
         End Do
      End Do
*
         End If   ! Seg1
*
         If (Seg2) Then
*
*-----Second quarter transformation
*
      Do iCntr2 = 1, nCntr2
         Do iZeta = 1, lZeta
            iPrm2 = (Indij(iZeta)-1)/nPrm1+1
            iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
*-----------Check for zero due to segmented basis
            If (Abs(Coeff2(iPrm2,iCntr2)).gt.Zero) Then
                  Do iVec = iiVec, iiVec+mVec-1
                     ijVec = mVec*(iPrm1-1) + (iVec-iiVec+1)
                     A3(iVec,iZeta) = A3(iVec,iZeta) +
     &                 Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec)
                  End Do
            End If
         End Do
      End Do
*
         Else
*
*-----Second quarter transformation
*
      Do iCntr2 = 1, nCntr2
         Do iZeta = 1, lZeta
            iPrm2 = (Indij(iZeta)-1)/nPrm1+1
            iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
               Do iVec = iiVec, iiVec+mVec-1
                  ijVec = mVec*(iPrm1-1) + (iVec-iiVec+1)
                  A3(iVec,iZeta) = A3(iVec,iZeta) +
     &              Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec)
               End Do
         End Do
      End Do
*
      End If
*
*-----End of loop sectioning
*
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(mZeta)
      End
