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

!#define _DEBUGPRINT_
subroutine MinDns(Dens,mBT,NumD,XCff,ltXCff,nD)
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

implicit none
integer mBT, nD, NumD, ltXCff
real*8 XCff(ltXCff,nD)
real*8, target :: Dens(mBT,nD,NumD)
! Define local variables
integer iC, iCol, iD, iM, iMat, iR, iRow, iStart, jCol, jRow
#ifdef _DEBUGPRINT_
integer i
#endif
real*8 XC
real*8 BVec(MxIter,2)
real*8, dimension(:,:,:), allocatable :: AMat
real*8, dimension(:,:), allocatable, target :: DRow, DCol
real*8, dimension(:,:), pointer :: pDR, pDC
real*8, external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call mma_allocate(DRow,nBT,nD,Label='DRow')
call mma_allocate(DCol,nBT,nD,Label='DCol')
call mma_allocate(AMat,MxIter,MxIter,2,label='AMat')
XCff(:,:) = Zero
AMat(:,:,:) = Zero
BVec(:,:) = Zero

! Allow a maximum of 10 densities in the minimization process to
! improve the numerical accuracy. This will also reduce the I/O.

!iStart = iter-iDMin
iStart = max(1,iter-9)
#ifdef _DEBUGPRINT_
write(6,*) 'iter,iStart=',iter,iStart
#endif

jRow = 0
do iRow=iStart,iter-1
  jRow = jRow+1

  iR = MapDns(iRow)
  if (iR < 0) then
    call RWDTG(-iR,DRow,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    pDR => DRow
  else
    pDR => Dens(1:mBT,1:nD,iR)
  end if

# ifdef _DEBUGPRINT_
  call NrmClc(pDR,nBT*nD,'MinDns','pDR')
# endif
  do iD=1,nD
    AMat(jRow,jRow,iD) = DDot_(nBT,pDR(:,iD),1,pDR(:,iD),1)
    BVec(jRow,iD) = DDot_(nBT,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
  end do ! iD

  jCol = 0
  do iCol=iStart,iRow-1
    jCol = jCol+1

    iC = MapDns(iCol)
    if (iC < 0) then
      call RWDTG(-iC,DCol,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
      pDC => DCol
    else
      pDC => Dens(1:mBT,1:nD,iC)
    end if

    do iD=1,nD
      AMat(jRow,jCol,iD) = DDot_(nBT,pDR(:,iD),1,pDC(:,iD),1)
      AMat(jCol,jRow,iD) = AMat(jRow,jCol,iD)
    end do ! iD

    nullify(pDC)

  end do ! iCol
  nullify(pDR)
end do ! iRow

do iD=1,nD

  ! Remove linear dependences from A-matrix
  call RmLDep(AMat(1,1,iD),MxIter,iter-iStart)

  ! Get minimization coefficients
  call DGEMM_('N','N',iter-iStart,1,iter-iStart, &
              One,AMat(1,1,iD),MxIter, &
              BVec(1,iD),iter-iStart, &
              Zero,XCff(iStart,iD),iter-iStart)
# ifdef _DEBUGPRINT_
  write(6,*) ' Coefficients minimizing density difference:'
  write(6,'(5f16.8)') (XCff(i,iD),i=1,iter-1)
  write(6,*)
# endif

end do ! iD

! Construct minimized density
do iMat=iter-1,iStart,-1

  iM = MapDns(iMat)
  if (iM < 0) then
    call RWDTG(-iM,DRow,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    pDR => DRow
  else
    pDR => Dens(1:mBT,1:nD,iM)
  end if

  do iD=1,nD
    XC = -XCff(iMat,iD)
    call daxpy_(nBT,XC,pDR(:,iD),1,Dens(1,iD,iPsLst),1)
  end do ! iD

end do ! iMat

call mma_deallocate(AMat)
call mma_deallocate(DCol)
call mma_deallocate(DRow)

end subroutine MinDns
