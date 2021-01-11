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
      Subroutine C1DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,
     &                  nFix,iP,MinWdw)
************************************************************************
*                                                                      *
*         References:                                                  *
*           C1-DIIS: P. Csaszar and P. Pulay, J. Mol. Struc.           *
*                    114, 31-34 (1984).                                *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December 1994                                            *
************************************************************************
      use Slapaf_Parameters, only: iOptC
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 q(nInter,nIter+1), dq(nInter,nIter),
     &       H(nInter,nInter), g(nInter,nIter+1),
     &       error(nInter,nIter), B((nIter+1)*(nIter+1)), RHS(nIter+1)
      Integer   iP(nIter), iRc
      Real*8, Allocatable:: A(:)
*
*     Statement function
*
      ij(i,j,lda)= (j-1)*lda + i
*
      iRout = 114
      iPrint = nPrint(iRout)
*
      Call mma_allocate(A,nInter**2,Label='A')
      call dcopy_(nInter**2,H,1,A,1)
      iRc = 0
      call dpotrf_('U',nInter,A,nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'C1DIIS(DPOTRF): iRC=',iRC
         Call Abend()
      End If
*
*-----Compute the new set of error vectors
*
      call dcopy_(nInter*nIter,g,1,Error,1)
      iRc = 0
      Call DPOTRS('U',nInter,nIter,A,nInter,Error,nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'C1DIIS(DPOTRS): iRC=',iRC
         Call Abend()
      End If
      If (iPrint.ge.99) Call RecPrt(' Error vectors',' ',
     &    error,nInter,nIter)
*
*-----Set up small system of linear equations
*     If more error vectors than degrees of freedom
*     exclude those with large error.
*
      Do 400 i = 1, nIter
         iP(i) = i
 400  Continue
*
*
*-----Bubble sort index array with respect to the magnitude of the
*     error vector.
*
      Do 405 i = 1, nIter-1
         If (iAnd(iOptC,16).eq.16) Then
            Err1 = DDot_(nInter,Error(1,iP(i)),1,
     &                         Error(1,iP(i)),1)
         Else If (iAnd(iOptC,32).eq.32) Then
            Err1 = DDot_(nInter,Error(1,iP(i)),1,
     &                             g(1,iP(i)),1)
         Else If (iAnd(iOptC,64).eq.64) Then
            Err1 = DDot_(nInter,    g(1,iP(i)),1,
     &                             g(1,iP(i)),1)
         Else
            Err1 = Zero
            Call WarningMessage(2,' Illegal iOptC setting!')
            Call Abend()
         End If
         ii = i
         Do 406 j = i+1, nIter
            If (iAnd(iOptC,16).eq.16) Then
               Err2 = DDot_(nInter,Error(1,iP(j)),1,
     &                            Error(1,iP(j)),1)
            Else If (iAnd(iOptC,32).eq.32) Then
               Err2 = DDot_(nInter,Error(1,iP(j)),1,
     &                                g(1,iP(j)),1)
            Else If (iAnd(iOptC,64).eq.64) Then
               Err2 = DDot_(nInter,    g(1,iP(j)),1,
     &                                g(1,iP(j)),1)
            Else
               Err2=Zero
               Call WarningMessage(2,' Illegal iOptC setting!')
               Call Abend()
            End If
            If (Err2.gt.Err1) Then
               ii = j
               Err1 = Err2
            End If
 406     Continue
         If (ii.ne.i) Then
            iSave = iP(i)
            iP(i) = iP(ii)
            iP(ii) = iSave
         End If
 405  Continue
      If (iPrint.ge.99) Write (6,*) ' iP=',iP
*
      MaxWdw=Max(2,(nInter-nFix)/2)
      mIter=Min(nIter,Min(MinWdw,MaxWdw))
      iOff = Max(0,nIter-mIter)
      B(ij(mIter+1,mIter+1,mIter+1)) = Zero
      RHS(mIter+1)       = -One
      Do 100 i = 1, mIter
         Do 110 j = 1, i-1
            If (iAnd(iOptC,16).eq.16) Then
               B(ij(i,j,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,
     &                                          error(1,iP(j+iOff)),1)
            Else If (iAnd(iOptC,32).eq.32) Then
               B(ij(i,j,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,
     &                                              g(1,iP(j+iOff)),1)
            Else If (iAnd(iOptC,64).eq.64) Then
               B(ij(i,j,mIter+1)) = DDot_(nInter,    g(1,iP(i+iOff)),1,
     &                                              g(1,iP(j+iOff)),1)
            Else
               Call WarningMessage(2,' Illegal iOptC setting!')
               Call Abend()
            End If
            B(ij(j,i,mIter+1)) = B(ij(i,j,mIter+1))
 110     Continue
         If (iAnd(iOptC,16).eq.16) Then
            B(ij(i,i,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,
     &                                       error(1,iP(i+iOff)),1)
         Else If (iAnd(iOptC,32).eq.32) Then
            B(ij(i,i,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,
     &                                           g(1,iP(i+iOff)),1)
         Else If (iAnd(iOptC,64).eq.64) Then
            B(ij(i,i,mIter+1)) = DDot_(nInter,    g(1,iP(i+iOff)),1,
     &                                           g(1,iP(i+iOff)),1)
         Else
            Call WarningMessage(2,' Illegal iOptC setting!')
            Call Abend()
         End If
         B(ij(i,mIter+1,mIter+1)) = -One
         B(ij(mIter+1,i,mIter+1)) = -One
         RHS(i)       = Zero
 100  Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' The B Matrix',' ',B,mIter+1,mIter+1)
         Call RecPrt(' The RHS',' ',RHS,1,mIter+1)
      End If
*
*-----Solve linear equation system
*
      Call Gauss(mIter+1,mIter+1,B,RHS,RHS)
      If (iPrint.ge.99) Call RecPrt(' The solution vector',
     &   ' ',RHS,1,mIter+1)
*
*-----Compute the interpolated parameter vector and
*     the interpolated gradient vector.
*
      call dcopy_(nInter,[Zero],0,q(1,nIter+1),1)
      call dcopy_(nInter,[Zero],0,g(1,nIter+1),1)
      Do 200 jIter = 1, mIter
         Call DaXpY_(nInter,RHS(jIter),q(1,iP(jIter+iOff)),1,
     &              q(1,nIter+1),1)
         Call DaXpY_(nInter,RHS(jIter),g(1,iP(jIter+iOff)),1,
     &              g(1,nIter+1),1)
 200  Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' The ipv',' ',q(1,nIter+1),1,nInter)
         Call RecPrt(' The igv',' ',g(1,nIter+1),1,nInter)
      End If
*
*-----Compute a new independent geometry by relaxation of
*     the interpolated gradient vector.
*
      call dcopy_(nInter,g(1,nIter+1),1,dq(1,nIter),1)
      Call DPOTRS('U',nInter,1,A,nInter,dq(1,nIter),nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'C1DIIS(DPOTRS): iRC=',iRC
         Call Abend()
      End If
      If (iPrint.ge.99) Call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
*
*     The shift is relative to the interpolated parameter
*     vector and we have to change it so that it is relative to the
*     actual parameter vector.
*
      Do 600 iInter = 1, nInter
         dq(iInter,nIter) =  dq(iInter,nIter) + q(iInter,nIter+1) -
     &                                          q(iInter,nIter  )
 600  Continue
      If (iPrint.ge.99) Call RecPrt(' dq(corr.)',' ',
     &   dq(1,nIter),1,nInter)
*
      Call mma_deallocate(A)
      Return
      End
