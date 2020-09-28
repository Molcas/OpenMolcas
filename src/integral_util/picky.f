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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      SubRoutine Picky(Din,n,m,nijPrm,nijCmp,nDCR,nSt,nEnd,mSt,mEnd,
     &                 Dout)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Din(((n*m+1)*nijCmp)+nijPrm+1,nDCR),
     &       Dout((((nEnd-nSt+1)*(mEnd-mSt+1)+1)*nijCmp)+nijPrm+1, nDCR)
*
*-----Statement function
*
      ij1(i,j,k) = (k-1)*(n*m+1) + (j-1)*n + i
      ij2(i,j,k) = (k-1)*((nEnd-nSt+1)*(mEnd-mSt+1)+1) +
     &             (j-1)*(nEnd-nSt+1) + i
*
      iRout = 236
      iPrint = nPrint(iRout)
*     Call GetMem('Enter Picky','Check','Real',iDum,iDum)
*
      If (nSt.eq.1.and.nEnd.eq.n .and.
     &    mSt.eq.1.and.mEnd.eq.m) Then
*--------Copy the whole block
         call dcopy_((((n*m+1)*nijCmp)+nijPrm+1)*nDCR,Din,1,Dout,1)
      Else
         iIn = (n*m+1)*nijCmp
         iOut = (((nEnd-nSt+1)*(mEnd-mSt+1)+1)*nijCmp)
*--------Loop over desymmetrized density blocks
         Do 1 iDCR = 1, nDCR
*-----------Loop over angular combinations
            Do 10 ijCmp= 1, nijCmp
               jm = 0
*--------------Loop over subset of contracted basis
               Do 20 im = mSt, mEnd
                  jm = jm + 1
                  jn = 0
*-----------------Loop over subset of contracted basis
                  Do 30 in = nSt, nEnd
                     jn = jn + 1
                     Dout(ij2(jn,jm,ijCmp),iDCR) =
     &                   Din(ij1(in,im,ijCmp),iDCR)
 30               Continue
 20            Continue
*--------------Move the largest density matrix element for
*              this angular combination
c              jIn = nijCmp*(n*m+1)
c              jOut= nijCmp*((nEnd-nSt+1)*(mEnd-mSt+1)+1)
               jIn = ij1(n,m,ijCmp)+1
               jOut= ij2(nEnd-nSt+1,mEnd-mSt+1,ijCmp)+1
               Dout(jOut,iDCR) = Din(jIn,iDCR)
 10         Continue
*-----------Move the largest density matrix element for
*           each pair plus the overall largest element
            call dcopy_(nijPrm+1,Din(iIn+1,iDCR),1,
     &                         Dout(iOut+1,iDCR),1)
 1       Continue
      End If
*
*     Call GetMem('Exit Picky','Check','Real',iDum,iDum)
      Return
      End
