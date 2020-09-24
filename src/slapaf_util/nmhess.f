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
      SubRoutine NmHess(dq,nInter,g,nIter,Hess,Delta,q,FEq,Cubic,
     &                  DipM,dDipM)
************************************************************************
*                                                                      *
* Object: to numerically evaluate the molecular Hessian.               *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             May '92                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 dq(nInter,nIter), g(nInter,nIter),Hess(nInter,nInter),
     &       q(nInter,nIter+1), FEq(nInter,nInter,nInter),
     &       DipM(3,nIter), dDipM(3,nInter)
      Logical Cubic
#include "print.fh"
#include "real.fh"
*
      iRout = 181
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt('NmHess:  g',' ',g,nInter,nIter)
         Call RecPrt('NmHess:  q',' ',q,nInter,nIter)
         Call RecPrt('NmHess: dq',' ',dq,nInter,nIter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Form the derivative of the dipole moment
*
      Do iInter = 1, nInter
         kIter  = 1 + (iInter-1)*2
         kIter1 = kIter + 1
         kIter2 = kIter + 2
C        Write (*,*) kIter1,kIter2
         dDipM(1,iInter) = (DipM(1,kIter1)-DipM(1,kIter2))
     &                   / (Two*Delta)
         dDipM(2,iInter) = (DipM(2,kIter1)-DipM(2,kIter2))
     &                   / (Two*Delta)
         dDipM(3,iInter) = (DipM(3,kIter1)-DipM(3,kIter2))
     &                   / (Two*Delta)
      End Do
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('DipM',' ',DipM,3,nIter)
      Call RecPrt('dDipM',' ',dDipM,3,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Evaluate the Hessian
*
      Do iInter = 1, nInter
         Do jInter = 1, nInter
            kIter = 1 + (jInter-1)*2
            kIter1 = kIter + 1
            kIter2 = kIter + 2
*           Write (*,*) iInter,jInter,kIter1,kIter2
*-----------Observe the sign convention due to the use of forces
*           rather than gradients!!!
            Hess(iInter,jInter) = -(g(iInter,kIter1) -
     &                              g(iInter,kIter2)) / (Two*Delta)
         End Do
      End Do
      If (iPrint.ge.99) Call RecPrt(' Numerical Hessian',' ',
     &   Hess,nInter,nInter)
*
*-----Symmetrize
*
      Do iInter = 1, nInter
         Do jInter = 1, iInter-1
            Hess(iInter,jInter) = (Hess(iInter,jInter) +
     &                             Hess(jInter,iInter) ) / Two
            Hess(jInter,iInter) =  Hess(iInter,jInter)
         End Do
      End Do
      If (iPrint.ge.99) Call RecPrt(' Symmetrized Hessian',' ',
     &   Hess,nInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the numerical cubic constants if data is available.
*
      If (Cubic) Then
*
*------- F(i,j,j)
*
         Do iInter = 1, nInter
            Do jInter = 1, nInter
               kIter = 1 + (jInter-1)*2
               kIter1 = kIter + 1
               kIter2 = kIter + 2
*------------- Observe the sign convention due to the use of forces
*              rather than gradients!!!
               FEq(iInter,jInter,jInter) = -(g(iInter,kIter1)
     &                                   +   g(iInter,kIter2))
     &                                   /  (Delta**2)
            End Do
         End Do
*
*------- F(i,j,k); j>k
*
         Do iInter = 1, nInter
            iCount=0
            Do jInter = 1, nInter
               Do kInter = 1, jInter-1
                  iCount=iCount+1
                  kIter = 1 + 2*nInter + (iCount-1)*4
                  kIter1= kIter + 1
                  kIter2= kIter + 2
                  kIter3= kIter + 3
                  kIter4= kIter + 4
*---------------- Observe the sign convention due to the use of forces
*                 rather than gradients!!!
                  FEq(iInter,jInter,kInter) = -(g(iInter,kIter1)
     &                                      -   g(iInter,kIter2)
     &                                      -   g(iInter,kIter3)
     &                                      +   g(iInter,kIter4))
     &                                      /   (Two*Delta)**2
               End Do
            End Do
         End Do
*
*------- Symmetrize
*
         Do iInter = 1, nInter
            Do jInter = 1, iInter
               Do kInter = 1, jInter
                  FEq(iInter,jInter,kInter) = (
     &               FEq(iInter,jInter,kInter) +
     &               FEq(iInter,kInter,jInter) +
     &               FEq(jInter,iInter,kInter) +
     &               FEq(jInter,kInter,iInter) +
     &               FEq(kInter,jInter,iInter) +
     &               FEq(kInter,iInter,jInter) ) / Six
                  FEq(iInter,jInter,kInter) =
     &               FEq(iInter,jInter,kInter)
                  FEq(iInter,kInter,jInter) =
     &               FEq(iInter,jInter,kInter)
                  FEq(jInter,iInter,kInter) =
     &               FEq(iInter,jInter,kInter)
                  FEq(jInter,kInter,iInter) =
     &               FEq(iInter,jInter,kInter)
                  FEq(kInter,iInter,jInter) =
     &               FEq(iInter,jInter,kInter)
               End Do
            End Do
         End Do
*
      End If
*
      Return
      End
