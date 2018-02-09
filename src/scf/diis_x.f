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
* Copyright (C) 1994,2017, Roland Lindh                                *
*               1995, Per-Olof Widmark                                 *
*               1995, Markus P. Fuelscher                              *
*               1995, Piotr Borowski                                   *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine DIIS_x(nD,CInter,nCI,QNRStp,HDiag,mOV,Ind)
************************************************************************
*                                                                      *
*     purpose: Accelerate convergence using DIIS method                *
*                                                                      *
*     input:                                                           *
*                                                                      *
*     output:                                                          *
*       CInter  : Interpolation coefficients of length nCI             *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: GrdClc, Gauss, OptClc                                  *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     modified by:                                                     *
*     P.O. Widmark, M.P. Fuelscher, P. Borowski & M.Schuetz            *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*     Derived from code for c1- and c2-DIIS as implemented by          *
*     R. Lindh in Slapaf and SCF in 1994.                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "llists.fh"
*
      Real*8 CInter(nCI,nD), HDiag(mOV,nD)
      Real*8, Dimension(:,:), Allocatable:: EVector, Bij
      Real*8, Dimension(:), Allocatable:: EValue, Err1, Err2, Scratch
*
*---- Define local variables
      Integer Ind(MxOptm)
      Real*8 GDiis(MxOptm + 1),BijTri(MxOptm*(MxOptm + 1)/2)
      Logical QNRstp
*define _NEW_CODE_
#ifdef _NEW_CODE_
      Logical Ignore
#endif
      Character*80 Text,Fmt
*
*---- Statement function for triangular index
      iTri(i,j) = i*(i-1)/2 + j
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*
#ifdef _NEW_CODE_
      Do i = kOptim, 1, -1
         Ind(i)=0
*
         tmp0= 0.0D0
         Do j = 1, iter
*
            Ignore=.False.
            Do k = kOptim, i+1, -1
               Ignore = Ignore .or. Ind(k).eq.j
            End Do
            If (Ignore) Cycle
*
            tmp = 0.0D0
            Do iD = 1, nD
               tmp = tmp + Elst(j,iD)
            End Do
*
            If (tmp.lt.tmp0) Then
               tmp0=tmp
               Ind(i)=j
            End If
*
         End Do
      End Do
      Write (6,*) (Ind(i),i=1,kOptim)
#else
      Do i = 1, kOptim
         Ind(i) = iter-kOptim+i
      End Do
#endif
*
*-----The following piece of code computes the DIIS coeffs
*     with error vectors chosen as the grds (DIIS only)
*     or as delta=-Hinv*grd (QNR/DIIS)
*
*-----in case of QNR step: grd already calculated in LinSer
*     or above...
*     in case of DIIS only: grd already calculated above
*
*     Allocate memory for error vectors (gradient or delta)
*
      Call mma_allocate(Err1,nOV*nD,Label='Err1')
      Call mma_allocate(Err2,nOV*nD,Label='Err2')
      nBij=kOptim+1
      Call mma_allocate(Bij,nBij,nBij)
      Call FZero(Bij,nBij**2)
*
*---- Compute norms, <e_i|e_j>
*
      Do i=1,kOptim
         Call ErrV(nOV*nD,Ind(i),QNRstp,Err1,HDiag)
*define _DEBUG_
#ifdef _DEBUG_
         Call NrmClc(Err1,nOV*nD,'Diis  ','Err1  ')
#endif
         Do j=1,i-1

            Call ErrV(nOV*nD,Ind(j),QNRstp,Err2,HDiag)
#ifdef _DEBUG_
            Call NrmClc(Err2,nOV*nD,'Diis  ','Err2  ')
#endif
            Bij(i,j) = DBLE(nD)*DDot_(nOV*nD,Err1,1,Err2,1)
            Bij(j,i) = Bij(i,j)
         End Do
         Bij(i,i) = DBLE(nD)*DDot_(nOV*nD,Err1,1,Err1,1)
      End Do
*
*---- Deallocate memory for error vectors & gradient
      Call mma_deallocate(Err2)
      Call mma_deallocate(Err1)
*
*---- Should we perform emergency convergence?
*
      Bsmall=Bij(1,1)
      Do i=2,kOptim
         Bsmall=Min(Bsmall,Bij(i,i))
      End Do
*
      If (Bsmall.lt.1.0d-30) Then
*        EmConv=.true.
*        Write (6,*) 'Bsmall.lt.1.0D-30'
*        Write (6,*) 'Bsmall=',BSmall
      End If
*
#ifdef _DEBUG_
      Write(6,*)'   Calculation of the norms in Diis :'
      Fmt  = '(6f16.8)'
      Text = 'B-matrix squared in Diis :'
      Call RecPrt(Text,Fmt,Bij,nBij,nBij)
      Write(6,*)
      Write(6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Here, the DIIS coeffs for DIIS only or QNR/DIIS are
*     computed, either with C1DIIS or C2DIIS
*     (-> stored in vector CInter)
*                                                                      *
************************************************************************
*                                                                      *
      If (.not.c1Diis) Then
*                                                                      *
*-------   C2DIIS case                                                 *
*                                                                      *
*         References:                                                  *
*         H. Sellers, Int. J. Quantum Chem. 45, 31-41(1993).           *
*         DOI: 10.1002/qua.560450106                                   *
*                                                                      *
************************************************************************
*                                                                      *
         If (QNRStp) Then
            AccCon = 'QNRc2DIIS'
         Else
            AccCon = 'c2DIIS   '
         End If
*
*------- Form a unit eigenvector matrix
         Call mma_allocate(EVector,kOptim,kOptim)
         Call mma_allocate(EValue,kOptim)
*
         call dcopy_(kOptim**2,Zero,0,EVector,       1)
         call dcopy_(kOptim,   One, 0,EVector,kOptim+1)
*
*------- Form a triangular B-matrix
*
         ij = 1
         Do i = 1, kOptim
            call dcopy_(i,Bij(i,1),nBij,BijTri(ij),1)
            ij = ij + i
         End Do
*
#ifdef _DEBUG_
         Fmt  = '(5g25.15)'
         Text = 'B-matrix before Jacobi :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenVectors before Jacobi :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*
*------- Diagonalize B-matrix
*
         EMax=0.0D0
         Do i = 1, kOptim*(kOptim+1)/2
            EMax=Max(EMax,Abs(BijTri(i)))
         End Do
         Do i = 1, kOptim*(kOptim+1)/2
            If (Abs(BijTri(i)).lt.EMax*1.0D-14) BijTri(i)=0.0D0
         End Do
*
         Call mma_allocate(Scratch,kOptim**2,Label='Scratch')
*
         Dummy=0.0D0
         iDum=0
         Call Diag_Driver('V','A','L',kOptim,BijTri,
     &                    Scratch,kOptim,Dummy,Dummy,iDum,iDum,
     &                    EValue,EVector,kOptim,1,0,'J',
     &                    nFound,iErr)
*
         Call mma_deallocate(Scratch)
         Call dCopy_(kOptim*(kOptim+1)/2,Zero,0,BijTri,1)
*
         iDiag = 0
         Do i = 1,kOptim
            iDiag = iDiag + i
            BijTri(iDiag) = EValue(i)
         End Do
*
#ifdef _DEBUG_
         Fmt  = '(5g25.15)'
         Text = 'B-matrix after Jacobi :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenValues :'
         Call RecPrt(Text,Fmt,EValue,1,kOptim)
         Text = 'EigenVectors :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*
*------  Renormalize the eigenvectors and eigenvalues to the
*        C1-DIIS format
*
#ifdef _DEBUG_
         Write(6,*)' Normalization constants :'
#endif
*
         Do iVec = 1, kOptim
            Alpha = Zero
            Do i = 1, kOptim
               Alpha = Alpha + EVector(i,iVec)
            End Do
*
#ifdef _DEBUG_
            Fmt = '(A7,i2,A4,f16.8)'
            Write(6,Fmt)' Alpha(',iVec,') = ',Alpha
#endif
*
            Alpha = One/Alpha
            Call DScal_(kOptim,Alpha,EVector(1,iVec),1)
            iPos = iTri(iVec,iVec)
            BijTri(iPos) = BijTri(iPos) * Alpha**2
            EValue(iVec) = EValue(iVec)   * Alpha**2
         End Do
*
#ifdef _DEBUG_
         Fmt  = '(6e16.8)'
         Text = 'B-matrix after scaling :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenValues after scaling :'
         Call RecPrt(Text,Fmt,EValue,1,kOptim)
         Text = 'EigenVectors after scaling :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*------  Set up thresholds
         Thrld  = 1.0D-15
         ThrCff = 1.0d+00*Ten
*
*------  Select a vector.
         ee1   = 1.0D+72
         cDotV = 1.0D+72
         t1    = 0.0D0
         ipBst =-99999999
         Do iVec = 1, kOptim
*
*           Pick up eigenvalue (ee2) and the norm of the
*           eigenvector (c2).
*
            ee2 = BijTri(iTri(iVec,iVec))
            c2 = DDot_(kOptim,EVector(1,iVec),1,EVector(1,iVec),1)
            t2 = Abs(EVector(kOptim,iVec))/Sqrt(c2)
*
*---------  Reject if <e|e> is too low (round-off),
*           analys further.
*
            If (ee2.lt.Thrld) Then
#ifdef _DEBUG_
               Fmt  = '(A,i2,5x,g12.6)'
               Text = '<e|e> is low,         iVec, <e|e> = '
               Write(6,Fmt)Text(1:36),iVec,ee2
#endif
*
*------------  Reject if coefficients are too large (linear dep.).
*
               If (Sqrt(c2).gt.ThrCff) Then
#ifdef _DEBUG_
                  Fmt  = '(A,i2,5x,g12.6)'
                  Text = 'c**2 is too large,     iVec, c**2 = '
                  Write(6,Fmt)Text(1:36),iVec,c2
#endif
                  Go To 520
               End If
            End If
*
*---------  Reject if coefficients are too large (linear dep.).
*
            If (Sqrt(c2).gt.ThrCff*Two) Then
#ifdef _DEBUG_
               Fmt  = '(A,i2,5x,g12.6)'
               Text = 'c**2 is too large,     iVec, c**2 = '
               Write(6,Fmt)Text(1:36),iVec,c2
#endif
               Go To 520
            End If
*
*-----------Keep the best candidate
*
            If (ee2*Five.lt.ee1) Then
*-----------   New vector much lower eigenvalue.
               ee1   = ee2
               ipBst = iVec
               cDotV = c2
               t2 = t1
            Else If (ee2.le.ee1*Three) Then
*-----------   New vector is close to the old vector.
*              Selection based on relative weight of the last
*              density.
               If (t2.gt.t1*4.0d0) Then
*--------------   New vector much better relative weight.
                  ee1   = ee2
                  ipBst = iVec
                  cDotV = c2
                  t2 = t1
               Else If (t2*1.2d0.lt.t1) Then
*--------------   Vectors are close in relative weight too!
*                 Select on eigenvalue only
                  If (ee2.lt.ee1) Then
                     ee1   = ee2
                     ipBst = iVec
                     cDotV = c2
                     t2 = t1
                  End If
               End If
            End If
*
 520        Continue
*
         End Do
*
         If (ipBst.lt.1 .or. ipBst.gt.kOptim) Then
            Write(6,*)' No proper solution found in C2-DIIS !'
            Fmt  = '(6e16.8)'
            Text = 'EigenValues :'
            Call RecPrt(Text,Fmt,EValue,1,kOptim)
            Text = 'EigenVectors :'
            Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
            Call ErrTra
            Call Quit_OnConvError()
         End If
         call dcopy_(kOptim,EVector(1,ipBst),1,CInter(1,1),1)
*
#ifdef _DEBUG_
         Write(6,*)
         Write(6,*)' Selected root :',ipBst
         Write(6,'(A,f16.8)')'  c**2 =         ',cDotV
#endif
*
         Call mma_deallocate(EValue)
         Call mma_deallocate(EVector)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
*------- C1DIIS                                                        *
*                                                                      *
*         References:                                                  *
*         P. Csaszar and P. Pulay, J. Mol. Struc., 114, 31-34 (1984).  *
*         doi:10.1016/S0022-2860(84)87198-7                            *
*                                                                      *
************************************************************************
*                                                                      *
         If (QNRStp) Then
            AccCon = 'QNRc1DIIS'
         Else
            AccCon = 'c1DIIS   '
         End If
*
*        Set up the missing part of the matrix in Eq. (5) and the
*        vector on the RHS in the same equation. Note the sign change!
*
         Do i = 1, kOptim
            Bij(kOptim + 1,i) = - One ! note sign change
            Bij(i,kOptim + 1) = - One ! note sign change
            GDiis(i)          =   Zero
         End Do
         Bij(kOptim + 1,kOptim + 1) =   Zero
         GDiis(kOptim + 1)          = - One  ! note sign change
*
#ifdef _DEBUG_
         Write(6,*)' B matrix in DIIS_e:'
         Do i = 1, kOptim + 1
            Write(6,'(7f16.8)')(Bij(i,j),j = 1, kOptim + 1),
     &                          GDiis(i)
         End Do
#endif
*
*------- Add penalty function, new version
*
* ...    Pfact1: c_k tend to one
* ...    Pfact2: c_i (i!=k) tend to zero
* ...    Pfact3: c_k tend to zero
*
         Pfact1 = Max(0.0d0  * Bsmall, 0.0d0  )
         Pfact2 = Max(1.0d-2 * Bsmall, 1.0d-14)
         Pfact3 = Max(1.0d-2 * Bsmall, 1.0d-14)
*
         Do i=1,kOptim
            Bij(i,i)=Bij(i,i)+Pfact2
         End Do
         Bij(kOptim,kOptim) = Bij(kOptim,kOptim) + Pfact1
         Bij(kOptim,kOptim) = Bij(kOptim,kOptim) + Pfact3 -Pfact2
         GDiis(kOptim) = Pfact1
*
*------- Condition the B matrix
*
         B11 = Sqrt(Bij(1,1)*Bij(kOptim,kOptim))
         Do i = 1, kOptim
            Do j = 1, kOptim
               Bij(i,j) = Bij(i,j)/B11
            End Do
         End Do
*
*------- Solve for the coefficients, solve the equations.
*
         Call Gauss(kOptim + 1,nBij,Bij,CInter(1,1),GDiis)
*
*------- Normalize sum of interpolation coefficients
*
         Fact = Zero
         Do i = 1, kOptim
            Fact = Fact + CInter(i,1)
         End Do
*
         Fact = One/Fact
         Do i = 1, kOptim
            CInter(i,1) = Fact*CInter(i,1)
         End Do
*
*------- Make sure new density gets a weight
         Call C_Adjust(CInter(1,1),kOptim,0.05D0)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*
      Call mma_deallocate(Bij)
*
*     Temporary fix for UHF.
*
      If (nD.eq.2) Call DCopy_(nCI,CInter(1,1),1,CInter(1,2),1)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Fmt  = '(6e16.8)'
      Text = 'The solution vector :'
      Call RecPrt(Text,Fmt,CInter(1,1),1,kOptim)
#endif
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 6) = TimFld( 6) + (Cpu2 - Cpu1)
      Return
      End
