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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

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

use InfSCF, only: iDisk, iPsLst, Iter, MapDns, MxIter, nBT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: mBT, NumD, ltXCff, nD
real(kind=wp), target, intent(inout) :: Dens(mBT,nD,NumD)
real(kind=wp), intent(out) :: XCff(ltXCff,nD)
integer(kind=iwp) :: iC, iCol, iD, iM, iMat, iR, iRow, iStart, jCol, jRow
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: BVec(MxIter,2)
real(kind=wp), allocatable :: AMat(:,:,:)
real(kind=wp), allocatable, target :: DRow(:,:), DCol(:,:)
real(kind=wp), pointer :: pDR(:,:), pDC(:,:)
real(kind=wp), external :: DDot_

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
write(u6,*) 'iter,iStart=',iter,iStart
#endif

jRow = 0
do iRow=iStart,iter-1
  jRow = jRow+1

  iR = MapDns(iRow)
  if (iR < 0) then
    call RWDTG(-iR,DRow,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    pDR => DRow
  else
    pDR => Dens(:,:,iR)
  end if

# ifdef _DEBUGPRINT_
  call NrmClc(pDR,nBT*nD,'MinDns','pDR')
# endif
  do iD=1,nD
    AMat(jRow,jRow,iD) = DDot_(nBT,pDR(:,iD),1,pDR(:,iD),1)
    BVec(jRow,iD) = DDot_(nBT,pDR(:,iD),1,Dens(:,iD,iPsLst),1)
  end do ! iD

  jCol = 0
  do iCol=iStart,iRow-1
    jCol = jCol+1

    iC = MapDns(iCol)
    if (iC < 0) then
      call RWDTG(-iC,DCol,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
      pDC => DCol
    else
      pDC => Dens(:,:,iC)
    end if

    do iD=1,nD
      AMat(jRow,jCol,iD) = DDot_(nBT,pDR(:,iD),1,pDC(:,iD),1)
    end do ! iD
    AMat(jCol,jRow,1:nD) = AMat(jRow,jCol,1:nD)

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
  write(u6,*) ' Coefficients minimizing density difference:'
  write(u6,'(5f16.8)') (XCff(i,iD),i=1,iter-1)
  write(u6,*)
# endif

end do ! iD

! Construct minimized density
do iMat=iter-1,iStart,-1

  iM = MapDns(iMat)
  if (iM < 0) then
    call RWDTG(-iM,DRow,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    pDR => DRow
  else
    pDR => Dens(:,:,iM)
  end if

  do iD=1,nD
    Dens(:,iD,iPsLst) = Dens(:,iD,iPsLst)-XCff(iMat,iD)*pDR(1:nBT,iD)
  end do ! iD

end do ! iMat

call mma_deallocate(AMat)
call mma_deallocate(DCol)
call mma_deallocate(DRow)

end subroutine MinDns
