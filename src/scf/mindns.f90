!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine MinDns(Dens,mBT,NumD,XCff,ltXCff,nD)
!***********************************************************************
!                                                                      *
!     purpose: Compute minimized density difference                    *
!                                                                      *
!     input:                                                           *
!       Dens    : a few last density matrix differences (mBT,NumD)     *
!                                                                      *
!     output:                                                          *
!       XCff    : coefficients (ltXCff==iDMin)                         *
!                                                                      *
!***********************************************************************
      use InfSCF, only: iDisk, iPsLst, Iter, nBT, MapDns
      use MxDM, only: MxIter
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer mBT, nD, NumD, ltXCff
      Real*8 XCff(ltXCff,nD)
      Real*8, Target:: Dens(mBT,nD,NumD)
!
!---- Define local variables
      Integer iC, iCol, iD, iM, iMat, iR, iRow, iStart, jCol, jRow
#ifdef _DEBUGPRINT_
      Integer i
#endif
      Real*8 XC
      Real*8 BVec(MxIter,2)
      Real*8, Dimension(:,:,:), Allocatable :: AMat
      Real*8, Dimension(:,:), Allocatable, Target:: DRow, DCol
      Real*8, Dimension(:,:), Pointer:: pDR, pDC
      Real*8, External:: DDot_
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!
      Call mma_allocate(DRow,nBT,nD,Label='DRow')
      Call mma_allocate(DCol,nBT,nD,Label='DCol')
      Call mma_allocate(AMat,MxIter,MxIter,2,label='AMat')
      XCff(:,:)=Zero
      AMat(:,:,:)=Zero
      BVec(:,:)=Zero
!
!     Allow a maximum of 10 densities in the minimization process to
!     improve the numerical accuracy. This will also reduce the I/O.
!
!     iStart = iter - iDMin
      iStart = Max(1,iter - 9)
#ifdef _DEBUGPRINT_
      Write (6,*) 'iter,iStart=',iter,iStart
#endif
!
      jRow=0
      Do iRow = iStart, iter - 1
         jRow = jRow + 1
!
         iR = MapDns(iRow)
         If (iR.lt.0) Then
            Call RWDTG(-iR,DRow,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
            pDR => DRow
         Else
            pDR => Dens(1:mBT,1:nD,iR)
         End If
!
#ifdef _DEBUGPRINT_
         Call NrmClc(pDR,nBT*nD,'MinDns','pDR')
#endif
         Do iD = 1, nD
            AMat(jRow,jRow,iD) = DDot_(nBT,pDR(:,iD),1,pDR(:,iD),1)
            BVec(jRow,iD) = DDot_(nBT,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
         End Do ! iD
!
         jCol = 0
         Do iCol = iStart, iRow - 1
            jCol = jCol + 1
!
            iC = MapDns(iCol)
            If (iC.lt.0) Then
               Call RWDTG(-iC,DCol,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
               pDC => DCol
            Else
               pDC => Dens(1:mBT,1:nD,iC)
            End If
!
            Do iD = 1, nD
               AMat(jRow,jCol,iD) = DDot_(nBT,pDR(:,iD),1,pDC(:,iD),1)
               AMat(jCol,jRow,iD) = AMat(jRow,jCol,iD)
            End Do ! iD
!
            Nullify(pDC)
!
         End Do ! iCol
         Nullify(pDR)
      End Do ! iRow
!
      Do iD = 1, nD
!
!----    Remove linear dependences from A-matrix
         Call RmLDep(AMat(1,1,iD),MxIter,iter-iStart)
!
!----    Get minimization coefficients
         Call DGEMM_('N','N',iter-iStart,1,iter-iStart,    &
                     One,AMat(1,1,iD),MxIter,            &
                           BVec(1,iD),iter-iStart,         &
                     Zero,XCff(iStart,iD),iter-iStart)
#ifdef _DEBUGPRINT_
         Write(6,*)' Coefficients minimizing density difference:'
         Write(6,'(5f16.8)')(XCff(i,iD),i=1,iter-1)
         Write(6,*)
#endif
!
      End Do ! iD
!
!---- Construct minimized density
      Do iMat = iter - 1, iStart, -1
!
         iM = MapDns(iMat)
         If (iM.lt.0) Then
            Call RWDTG(-iM,DRow,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
            pDR => DRow
         Else
            pDR => Dens(1:mBT,1:nD,iM)
         End If

         Do iD = 1, nD
            XC = - XCff(iMat,iD)
            call daxpy_(nBT,XC,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
         End Do ! iD
!
      End Do ! iMat
!
      Call mma_deallocate(AMat)
      Call mma_deallocate(DCol)
      Call mma_deallocate(DRow)

      end subroutine MinDns


      SubRoutine RmLDep(AMat,lDm,lth)
!***********************************************************************
!                                                                      *
!     purpose: Remove linear dependencies                              *
!                                                                      *
!     input:                                                           *
!       AMat    : matrix with linear dependencies (lDm,lDm)            *
!                                                                      *
!     output:                                                          *
!       AMat    : inverted matrix without linear dependencies (lDm,lDm)*
!                                                                      *
!***********************************************************************
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer lDm, lth
      Real*8 AMat(lDm,lDm)
!
!     Local variables
!
      Real*8 Dummy
      Integer i, iDum, iErr, ij, lthS, nFound, lthT
#ifdef _DEBUGPRINT_
      Integer j
#endif
      Real*8, Dimension(:), Allocatable:: ATri, EVec, EVal, Scr
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!
      lthT = lth*(lth + 1)/2
      lthS = lth*lth
      Call mma_allocate(ATri,lthT,Label='ATri')
      Call mma_allocate(EVec,lthS,Label='EVec')
      Call mma_allocate(EVal,lth ,Label='EVal')
!
!---- Put a unit matrix into the eigenvectors work space
      call dcopy_(lthS,[Zero],0,EVec,      1)
      call dcopy_(lth, [One], 0,EVec,lth + 1)
!
!---- Copy trialangular part of AMat to work space
      ij = 1
      Do i = 1, lth
         call dcopy_(i,AMat(i,1),lDm,ATri(ij),1)
         ij = ij + i
      End Do
#ifdef _DEBUGPRINT_
      Write(6,*)' Squared A-matrix in RmLDep:'
      Do i = 1, lth
         Write(6,'(5(1x,es13.6))')(AMat(i,j),j=1,lth)
      End Do
      Write(6,*)' Triangular A-matrix:'
      Write(6,'(5(1x,es13.6))')(ATri(i),i=1,lthT)
#endif
!
!---- Diagonalize
      Call mma_allocate(Scr,lth**2,Label='Scr')
      Dummy=Zero
      iDum=0
      Call Diag_Driver('V','A','L',lth,ATri,Scr,lth,        &
                       Dummy,Dummy,iDum,iDum,EVal,EVec,     &
                       lth,1,0,'J',nFound,iErr)
      Call mma_deallocate(Scr)
#ifdef _DEBUGPRINT_
      Write(6,*)' Eigenvalues of A-matrix in RLnDep:'
      Write(6,'(5(1x,es13.6))')(EVal(i),i=1,lth)
#endif
!
!---- Form the inverse
      Call dCopy_(lDm*lth,[Zero],0,AMat,1)
      Do i = 1, lth
         If (EVal(i).gt.1.0d-12) Then
            AMat(i,i) = 1.0d+00/EVal(i)
         Else
!           Write(*,*)' Eigenvalue',i,' smaller then 1.0e-12'
            AMat(i,i) = Zero
         End If
      End Do
      Call mma_allocate(Scr,lthS,Label='Scr')
      Call DGEMM_('N','T',lth,lth,lth,      &
                  One,AMat,lDm,             &
                        EVec,lth,           &
                  Zero,Scr,lth)
      Call DGEMM_('N','N',lth,lth,lth,      &
                  One,EVec,lth,             &
                        Scr,lth,            &
                  Zero,AMat,lDm)
      Call mma_deallocate(Scr)
!
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      Call mma_deallocate(ATri)
!
      End SubRoutine RmLDep
