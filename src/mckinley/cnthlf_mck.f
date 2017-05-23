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
      Subroutine Cnthlf_mck(Coeff1,nCntr1,nPrm1,
     &                      Coeff2,nCntr2,nPrm2,
     &                      nZeta,lZeta,nVec,First,IncVec,A1,A2,A3,
     &                  Indij)
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
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
c#include "print.fh"
      Real*8 Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2),
     &       A1(lZeta,nVec), A2(nPrm2,IncVec*nCntr1),
     &       A3(nVec,nCntr1,nCntr2)
      Integer Indij(nZeta)
      Logical First, Seg1, Seg2
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
      If (First) call dcopy_(nVec*nCntr1*nCntr2,Zero,0,A3,1)
*
*-----Loop sectioning
*
      Do iiVec = 1, nVec, IncVec
         mVec = Min(IncVec,nVec-iiVec+1)
*--------Set intermediate matrix to zero
         call dcopy_(nPrm2*nCntr1*mVec,Zero,0,A2,1)
*
         If (Seg1) Then
*
*-----First quarter transformation
*
      i = 0
*
      Do iPrm1 = 1, nPrm1
         Do iCntr1 = 1, nCntr1
*-----------Check for zero due to segmented basis
            If (Abs(Coeff1(iPrm1,iCntr1)).gt.Zero) Then
               Do iPrm2 = 1, nPrm2
                  iZeta = (iPrm2-1)*nPrm1 + iPrm1
                  jZeta = Indij(iZeta)
*-----------------Skip due to screening
                  If (jZeta.gt.0) Then
                     Do iVec = iiVec, iiVec+mVec-1
                        ijVec = mVec*(iCntr1-1) + (iVec-iiVec+1)
                        A2(iPrm2,ijVec) = A2(iPrm2,ijVec) +
     &                    Coeff1(iPrm1,iCntr1)*A1(jZeta,iVec)
                     End Do
                  End If
               End Do
            End If
         End Do
      End Do
*
         Else
*
*-----First quarter transformation
*
      i = 0
      Do iPrm1 = 1, nPrm1
         Do iCntr1 = 1, nCntr1
            Do iPrm2 = 1, nPrm2
               iZeta = (iPrm2-1)*nPrm1 + iPrm1
               jZeta = Indij(iZeta)
*--------------Skip due to screening
               If (jZeta.gt.0) Then
                  Do iVec = iiVec, iiVec+mVec-1
                     ijVec = mVec*(iCntr1-1) + (iVec-iiVec+1)
                     A2(iPrm2,ijVec) = A2(iPrm2,ijVec) +
     &                 Coeff1(iPrm1,iCntr1)*A1(jZeta,iVec)
                  End Do
               End If
            End Do
         End Do
      End Do
*
         End If
*
         If (Seg2) Then
*
*-----Second quarter transformation
*
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
*-----------Check for zero due to segmented basis
            If (Abs(Coeff2(iPrm2,iCntr2)).gt.Zero) Then
               Do iCntr1 = 1, nCntr1
                  Do iVec = iiVec, iiVec+mVec-1
                     ijVec = mVec*(iCntr1-1) + (iVec-iiVec+1)
                     A3(iVec,iCntr1,iCntr2) = A3(iVec,iCntr1,iCntr2) +
     &                 Coeff2(iPrm2,iCntr2)*A2(iPrm2,ijVec)
                  End Do
               End Do
            End If
         End Do
      End Do
*
         Else
*
*-----Second quarter transformation
*
      Do iPrm2 = 1, nPrm2
         Do iCntr2 = 1, nCntr2
            Do iCntr1 = 1, nCntr1
               Do iVec = iiVec, iiVec+mVec-1
                  ijVec = mVec*(iCntr1-1) + (iVec-iiVec+1)
                  A3(iVec,iCntr1,iCntr2) = A3(iVec,iCntr1,iCntr2) +
     &              Coeff2(iPrm2,iCntr2)*A2(iPrm2,ijVec)
               End Do
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
      End
