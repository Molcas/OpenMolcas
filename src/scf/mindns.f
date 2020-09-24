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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine MinDns(Dens,mBT,NumD,XCff,ltXCff,nD)
************************************************************************
*                                                                      *
*     purpose: Compute minimized density difference                    *
*                                                                      *
*     input:                                                           *
*       Dens    : a few last density matrix differences (mBT,NumD)     *
*                                                                      *
*     output:                                                          *
*       XCff    : coefficients (ltXCff==iDMin)                         *
*                                                                      *
*     called from: DMat, OptClc                                        *
*                                                                      *
*     calls to: RWDTG, Gauss                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Integer ltXCff
*
      Real*8, Target:: Dens(mBT,nD,NumD)
      Real*8 XCff(ltXCff,nD)
*
*---- Define local variables
      Real*8 BVec(MxIter,2)
      Real*8, Dimension(:,:,:), Allocatable :: AMat
      Real*8, Dimension(:,:), Allocatable, Target:: DRow, DCol
      Real*8, Dimension(:,:), Pointer:: pDR, pDC
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call qEnter('MinDns')
      Write(6,*)' ***** SubRoutine MinDns *****'
#endif
*
      Call mma_allocate(DRow,nBT,nD,Label='DRow')
      Call mma_allocate(DCol,nBT,nD,Label='DCol')
      Call mma_allocate(AMat,MxIter,MxIter,2,label='AMat')
      Call FZero(XCff,ltXCff*nD)
      Call FZero(AMat,2*MxIter**2)
      Call FZero(BVec,2*MxIter   )
*
      iter_d=iter-iter0
*
*     Allow a maximum of 10 densities in the minimization process to
*     improve the numerical accuracy. This will also reduce the I/O.
*
*     iStart = iter_d - iDMin
      iStart = Max(1,iter_d - 9)
#ifdef _DEBUGPRINT_
      Write (6,*) 'iter_d,iStart=',iter_d,iStart
#endif
*
      jRow=0
      Do iRow = iStart, iter_d - 1
         jRow = jRow + 1
*
         iR = MapDns(iRow)
         If (iR.lt.0) Then
            Call RWDTG(-iR,DRow,nBT*nD,'R','DENS  ',iDisk,MxDDsk)
            pDR => DRow
         Else
            pDR => Dens(1:mBT,1:nD,iR)
         End If
*
#ifdef _DEBUGPRINT_
         Call NrmClc(pDR,nBT*nD,'MinDns','pDR')
#endif
         Do iD = 1, nD
            AMat(jRow,jRow,iD) = DDot_(nBT,pDR(:,iD),1,pDR(:,iD),1)
            BVec(jRow,iD) = DDot_(nBT,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
         End Do ! iD
*
         jCol = 0
         Do iCol = iStart, iRow - 1
            jCol = jCol + 1
*
            iC = MapDns(iCol)
            If (iC.lt.0) Then
               Call RWDTG(-iC,DCol,nBT*nD,'R','DENS  ',iDisk,MxDDsk)
               pDC => DCol
            Else
               pDC => Dens(1:mBT,1:nD,iC)
            End If
*
            Do iD = 1, nD
               AMat(jRow,jCol,iD) = DDot_(nBT,pDR(:,iD),1,pDC(:,iD),1)
               AMat(jCol,jRow,iD) = AMat(jRow,jCol,iD)
            End Do ! iD
*
            Nullify(pDC)
*
         End Do ! iCol
         Nullify(pDR)
      End Do ! iRow
*
      Do iD = 1, nD
*
*----    Remove linear dependences from A-matrix
         Call RmLDep(AMat(1,1,iD),MxIter,iter_d-iStart)
*
*----    Get minimization coefficients
         Call DGEMM_('N','N',
     &               iter_d-iStart,1,iter_d-iStart,
     &               1.0d0,AMat(1,1,iD),MxIter,
     &                     BVec(1,iD),iter_d-iStart,
     &               0.0d0,XCff(iStart,iD),iter_d-iStart)
#ifdef _DEBUGPRINT_
         Write(6,*)' Coefficients minimizing density difference:'
         Write(6,'(5f16.8)')(XCff(i,iD),i=1,iter_d-1)
         Write(6,*)
#endif
*
      End Do ! iD
*
*---- Construct minimized density
      Do iMat = iter_d - 1, iStart, -1
*
         iM = MapDns(iMat)
         If (iM.lt.0) Then
            Call RWDTG(-iM,DRow,nBT*nD,'R','DENS  ',iDisk,MxDDsk)
            pDR => DRow
         Else
            pDR => Dens(1:mBT,1:nD,iM)
         End If

         Do iD = 1, nD
            XC = - XCff(iMat,iD)
            call daxpy_(nBT,XC,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
         End Do ! iD
*
      End Do ! iMat
*
      Call mma_deallocate(AMat)
      Call mma_deallocate(DCol)
      Call mma_deallocate(DRow)
*
#ifdef _DEBUGPRINT_
      Call qExit('MinDns')
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
      SubRoutine RmLDep(AMat,lDm,lth)
************************************************************************
*                                                                      *
*     purpose: Remove linear dependencies                              *
*                                                                      *
*     input:                                                           *
*       AMat    : matrix with linear dependencies (lDm,lDm)            *
*                                                                      *
*     output:                                                          *
*       AMat    : inverted matrix without linear dependencies (lDm,lDm)*
*                                                                      *
*     called from: MinDns                                              *
*                                                                      *
*     calls to:                                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 AMat(lDm,lDm)
*
*     Local variables
*
      Real*8, Dimension(:), Allocatable:: ATri, EVec, EVal, Scr
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      Call qEnter('RmLDep')
      Write(6,*)' ***** SubRoutine RmLDep*****'
#endif
*
      lthT = lth*(lth + 1)/2
      lthS = lth*lth
      Call mma_allocate(ATri,lthT,Label='ATri')
      Call mma_allocate(EVec,lthS,Label='EVec')
      Call mma_allocate(EVal,lth ,Label='EVal')
*
*---- Put a unit matrix into the eigenvectors work space
      call dcopy_(lthS,[Zero],0,EVec,      1)
      call dcopy_(lth, [One], 0,EVec,lth + 1)
*
*---- Copy trialangular part of AMat to work space
      ij = 1
      Do i = 1, lth
         call dcopy_(i,AMat(i,1),lDm,ATri(ij),1)
         ij = ij + i
      End Do
#ifdef _DEBUGPRINT_
      Write(6,*)' Squared A-matrix in RmLDep:'
      Do i = 1, lth
         Write(6,'(5(2x,e12.6))')(AMat(i,j),j=1,lth)
      End Do
      Write(6,*)' Triangular A-matrix:'
      Write(6,'(5(2x,e12.6))')(ATri(i),i=1,lthT)
#endif
*
*---- Diagonalize
      Call mma_allocate(Scr,lth**2,Label='Scr')
      Dummy=0.0D0
      iDum=0
      Call Diag_Driver('V','A','L',lth,ATri,Scr,lth,
     &                 Dummy,Dummy,iDum,iDum,EVal,EVec,
     &                 lth,1,0,'J',nFound,iErr)
      Call mma_deallocate(Scr)
#ifdef _DEBUGPRINT_
      Write(6,*)' Eigenvalues of A-matrix in RLnDep:'
      Write(6,'(5(2x,e12.6))')(EVal(i),i=1,lth)
#endif
*
*---- Form the inverse
      Call dCopy_(lDm*lth,[Zero],0,AMat,1)
      Do i = 1, lth
         If (EVal(i).gt.1.0d-12) Then
            AMat(i,i) = 1.0d+00/EVal(i)
         Else
*           Write(*,*)' Eigenvalue',i,' smaller then 1.0e-12'
            AMat(i,i) = Zero
         End If
      End Do
      Call mma_allocate(Scr,lthS,Label='Scr')
      Call DGEMM_('N','T',
     &            lth,lth,lth,
     &            1.0d0,AMat,lDm,
     &                  EVec,lth,
     &            0.0d0,Scr,lth)
      Call DGEMM_('N','N',
     &            lth,lth,lth,
     &            1.0d0,EVec,lth,
     &                  Scr,lth,
     &            0.0d0,AMat,lDm)
      Call mma_deallocate(Scr)
*
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      Call mma_deallocate(ATri)
*
#ifdef _DEBUGPRINT_
      Call qExit('RmLDep')
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
