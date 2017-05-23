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
* Copyright (C) 1994,1995, Roland Lindh                                *
************************************************************************
      Subroutine C2DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,
     &                  Scrt1,nScrt1,nFix,iP,iOptC)
************************************************************************
*                                                                      *
*         References:                                                  *
*           C2-DIIS: Sellers, Int. J. Quantum Chem. 45, 31-41(1993).   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
*                                                                      *
*             Modified for anharmonic constants by R. Lindh, Oct. '95  *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 q(nInter,nIter+1), dq(nInter,nIter),
     &       H(nInter,nInter), g(nInter,nIter+1),
     &       error(nInter,nIter+1), B((nIter+1)*(nIter+1)),
     &       RHS(nIter+1), Scrt1(nScrt1)
      Integer   iP(nIter), iRc
      Logical Fail
*
*     Statement function
*
      ij(i,j,lda)= (j-1)*lda + i
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      Call QEnter('C2DIIS')
      iRout = 121
      iPrint = nPrint(iRout)
*
      Call FZero(Error,nInter*(nIter+1))
*
      Call Allocate_Work(ipA,nInter**2)
      call dcopy_(nInter**2,H,1,Work(ipA),1)
      iRc = 0
      call dpotrf_('U',nInter,Work(ipA),nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'C2DIIS(DPOTRF): iRC=',iRC
         Call Abend()
      End If

*
*-----Compute the new set of error vectors. Here we have two
*     possibilities depending on if we use a 2nd or 3rd order
*     update method.
*
*     Note!!!!
*
*     i)  We store the force in g Not the gradient
*
*     ii) e is the displacement which should be added to the
*         current position to get to the equilibrium geometry.
*
*        r  = r   + e      e   =r  - r
*         eq   i-1   i-1    i-1  eq   i-1
*
*     Compute:
*        r  = r + e        e =r  - r
*         eq   i   i        k  eq   k
*
      If (iAnd(iOptC,16).eq.16 .or.
     &    iAnd(iOptC,32).eq.32) Then
         Call ThrdO(nInter,g(1,nIter),Work(ipA),Error(1,nIter),Fail)
         If (Fail) Then
            Call WarningMessage(2,'C2Diis: ThrdO Failed!')
            Call Quit_OnConvError()
         End If
         Do iIter = 1, nIter-1
            Do iInter = 1, nInter
               Error(iInter,iIter)=Error(iInter,nIter)
     &                            +q(iInter,nIter)-q(iInter,iIter)
            End Do
         End Do
      End If
      If (iPrint.ge.99) Call RecPrt(' Error vectors',' ',
     &    error,nInter,nIter)
*
*-----Set up small system of linear equations
*     If more error vectors than degrees of freedom
*     exclude those with large error.
*
      Do i = 1, nIter
         iP(i) = i
      End Do
*
*-----Bubble sort index array with respect to the magnitude of
*
*     <g|g>
*     <g|dx>
*     <dx|g>
*     <dx|dx>
*
      Do i = Max(1,nIter-11), nIter-1
         If (iAnd(iOptC,16).eq.16) Then
            Err1 = DDot_(nInter,Error(1,iP(i)),1,Error(1,iP(i)),1)
         Else If (iAnd(iOptC,32).eq.32) Then
            Err1 = DDot_(nInter,g(1,iP(i)),1,Error(1,iP(i)),1)
         Else If (iAnd(iOptC,64).eq.64) Then
            Err1 = DDot_(nInter,g(1,iP(i)),1,g(1,iP(i)),1)
         Else
            Err1=Zero
            Call WarningMessage(2,' Illegal iOptC setting!')
            Call Quit_OnUserError()
         End If
         ii = i
         Do j = i+1, nIter
            If (iAnd(iOptC,16).eq.16) Then
               Err2 = DDot_(nInter,Error(1,iP(j)),1,Error(1,iP(j)),1)
            Else If (iAnd(iOptC,32).eq.32) Then
               Err2 = DDot_(nInter,g(1,iP(j)),1,Error(1,iP(j)),1)
            Else If (iAnd(iOptC,64).eq.64) Then
               Err2 = DDot_(nInter,g(1,iP(j)),1,g(1,iP(j)),1)
            Else
               Err2=Zero
               Call WarningMessage(2,' Illegal iOptC setting!')
               Call Quit_OnUserError()
            End If
            If (Err2.gt.Err1) Then
               ii = j
               Err1 = Err2
            End If
         End Do
         If (ii.ne.i) Then
            iSave = iP(i)
            iP(i) = iP(ii)
            iP(ii) = iSave
         End If
      End Do
      If (iPrint.ge.99) Write (6,*) ' iP=',iP
*
*     MaxWdw=Max(3,3*(nInter-nFix)/4)
      MaxWdw=Max(3,(nInter-nFix)/2)
      MinWdw=Min(5,MaxWdw)
      mIter=Min(nIter,MinWdw)
      iOff = Max(0,nIter-mIter)
      Thrhld = 0.1D-13
      ThrCff = (Two*Ten)**2
      ThrLdp = Ten**3
      Do i = 1, mIter
         Do j = 1, i
            If (iAnd(iOptC,16).eq.16) Then
               B(iTri(i,j)) = DDot_(nInter,Error(1,iP(i+iOff)),1,
     &                                    Error(1,iP(j+iOff)),1)
            Else If (iAnd(iOptC,32).eq.32) Then
               B(iTri(i,j)) = DDot_(nInter,Error(1,iP(i+iOff)),1,
     &                                        g(1,iP(j+iOff)),1)
            Else If (iAnd(iOptC,64).eq.64) Then
               B(iTri(i,j)) = DDot_(nInter,    g(1,iP(i+iOff)),1,
     &                                        g(1,iP(j+iOff)),1)
            Else
               Call WarningMessage(2,' Illegal iOptC setting!')
               Call Quit_OnUserError()
            End If
         End Do
      End Do
      If (iPrint.ge.99) Call TriPrt(' The B Matrix',' ',B,mIter)
      call dcopy_(mIter**2,Zero,0,Scrt1,1)
      call dcopy_(mIter,One,0,Scrt1,mIter+1)
      Call NIDiag_new(B,Scrt1,mIter,mIter,0)
      If (iPrint.ge.99) Then
         Call TriPrt(' The B Matrix after diagonalization','(9E10.2)',
     &                B,mIter)
         Call RecPrt(' Eigenvectors','(9E10.2)',Scrt1,mIter,mIter)
      End If
*
*-----Renormalize the eigenvectors and eigenvalues to the
*     C1-DIIS format.
*
      Do iVec = 1, mIter
         Alpha = Zero
         Do i = 1, mIter
            Alpha = Alpha + Scrt1(ij(i,iVec,mIter))
         End Do
         Alpha = One/Alpha
         Call DScal_(mIter,Alpha,Scrt1(ij(1,iVec,mIter)),1)
         B(iTri(iVec,iVec)) = B(iTri(iVec,iVec)) * Alpha**2
      End Do
      If (iPrint.ge.99) Then
         Write (6,*) ' After normalization to C1-DIIS format'
         Call TriPrt(' The B Matrix after diagonalization','(9E10.2)',
     &              B,mIter)
         Call RecPrt(' Eigenvectors',' ',Scrt1,mIter,mIter)
      End If
*
*-----------Select a vector.
*
      ee_old = 1.0D+72
      c2_old = 1.0D+72
      iVec_old = -99999999
      Do iVec = 1, mIter
         If (iPrint.ge.99) Write (6,*) ' Scanning vector',iVec
         ee_new = B(iTri(iVec,iVec))
         If (iPrint.ge.99) Write (6,*)
     &           ' ee_old, ee_new=',ee_old, ee_new
*
*--------Examine if <e|e> is too low (possible round-off) or
*        linear dependency.
*
         If (ee_new.lt.Thrhld) Then
            If (iPrint.ge.99) Write (6,*)
     &         ' <e|e> is low in DIIS, iVec,<e|e>=',iVec,ee_new
*
*--------Reject if coefficients are too large (linear dep.).
*
            c2_new = DDot_(mIter,Scrt1(ij(1,iVec,mIter)),1,
     &                          Scrt1(ij(1,iVec,mIter)),1)
            If (c2_new.gt.ThrCff) Then
               If (iPrint.ge.99) Write (6,*)
     &            ' c**2 is too large in DIIS, iVec,c**2=',iVec,c2_new
                Go To 10
             End If
         End If
*
*--------Reject if coefficients are by far too large (linear dep.).
*
         c2_new = DDot_(mIter,Scrt1(ij(1,iVec,mIter)),1,
     &                   Scrt1(ij(1,iVec,mIter)),1)
         If (c2_new.gt.ThrLdp) Then
            If (iPrint.ge.99) Write (6,*)
     &         ' c**2 is too large in DIIS, iVec,c**2=',iVec,c2_new
            Go To 10
         End If
*
*--------Keep the best candidate
*
         If (ee_new*Five.lt.ee_old) Then
*-----------New vector much lower eigenvalue.
            c2_old=c2_new
            ee_old=ee_new
            iVec_old=iVec
            If (iPrint.ge.99) Write (6,*)
     &          'New vector much lower eigenvalue',iVec_old
         Else If (ee_new.le.ee_old*Five) Then
*-----------New vector is close to the old vector.
*           Selection based on relative weight of the last
*           geometry.
            If (iPrint.ge.99) Write (6,*)
     &          'Eigenvalues are close',iVec_old,iVec
            t1=Abs(Scrt1(ij(mIter,iVec_old,mIter)))/Sqrt(c2_old)
            t2=Abs(Scrt1(ij(mIter,iVec,    mIter)))/Sqrt(c2_new)
            If (t2.gt.t1*1.2d0) Then
*--------------New vector much better relative weight.
               c2_old=c2_new
               ee_old=ee_new
               iVec_old=iVec
               If (iPrint.ge.99) Write (6,*)
     &            'New vector much better relative weight',iVec_old
            Else If (t2*1.2d0.lt.t1) Then
*--------------Vectors are close in relative weight too!
*              Select on eigenvalue only
               If (iPrint.ge.99) Write (6,*)
     &            'Relative weights are close',iVec_old,iVec
               If (ee_new.lt.ee_old) Then
                  c2_old=c2_new
                  ee_old=ee_new
                  iVec_old=iVec
                  If (iPrint.ge.99) Write (6,*)
     &               'New vector has lower eigenvalue',iVec_old
               End If
            End If
         End If
*
 10      Continue
      End Do
      If (iVec_old.lt.1 .or. iVec_old.gt.mIter) Then
         Call WarningMessage(2,
     &               ' No proper solution found in C2-DIIS!')
         Call Abend()
      End If
      call dcopy_(mIter,Scrt1(ij(1,iVec_old,mIter)),1,RHS,1)
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Selecting root',iVec_old
         Call RecPrt(' The solution vector',' ',RHS,1,mIter)
      End If
*
*-----Compute the interpolated parameter vector and
*     the interpolated gradient vector.
*
      call dcopy_(nInter,Zero,0,q(1,nIter+1),1)
      call dcopy_(nInter,Zero,0,g(1,nIter+1),1)
      call dcopy_(nInter,Zero,0,Scrt1,1)
      Do iIter = 1, mIter
*
         If (Abs(RHS(iIter)).lt.1.0D-12) Go To 11
*
         Call DaXpY_(nInter,RHS(iIter),Error(1,iP(iIter+iOff)),1,
     &              Scrt1,1)
*
*------- The interpolated parameter vector is computed as a
*        simple linear combination
*
         Call DaXpY_(nInter,RHS(iIter),q(1,iP(iIter+iOff)),1,
     &              q(1,nIter+1),1)
*
*------- The interpolated gradient vector (Stored as force)
*
*        Sum(i) c  g                 (Coeffs stored in RHS)
*                i  i
*
         Call DaXpY_(nInter,RHS(iIter),g(1,iP(iIter+iOff)),1,
     &              g(1,nIter+1),1)
*
 11      Continue
      End Do

      If (iPrint.ge.99) Then
         Call RecPrt(' The iev',' ',Scrt1,1,nInter)
         Call RecPrt(' The ipv',' ',q(1,nIter+1),1,nInter)
         Call RecPrt(' The igv',' ',g(1,nIter+1),1,nInter)
      End If
*
*-----Compute a new independent geometry by relaxation of
*     the interpolated gradient vector.
*
      call dcopy_(nInter,g(1,nIter+1),1,dq(1,nIter),1)
      iRc = 0
      Call DPOTRS('U',nInter,1,Work(ipA),nInter,dq(1,nIter),nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'C2DIIS(DPOTRS): iRC=',iRC
         Call Abend()
      End If
      If (iPrint.ge.99) Call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
*
*     The shift is relative to the interpolated parameter
*     vector and we have to change it so that it is relative to the
*     actual parameter vector.
*
      Do iInter = 1, nInter
         dq(iInter,nIter) =  dq(iInter,nIter) + q(iInter,nIter+1) -
     &                                          q(iInter,nIter  )
      End Do
      If (iPrint.ge.99) Call RecPrt(' dq(corr.)',' ',
     &   dq(1,nIter),1,nInter)
*
      Call Free_Work(ipA)
*
      Call QExit('C2DIIS')
      Return
      End
