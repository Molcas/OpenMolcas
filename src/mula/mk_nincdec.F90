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
! Copyright (C) 1995,1998, Niclas Forsberg                             *
!               1995,1998, Anders Bernhardsson                         *
!               2009, Giovanni Ghigo                                   *
!***********************************************************************

subroutine Mk_nIncDec(m_max,nOrd,msiz,mInc,mDec,mMat,Graph2,nOsc)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: m_max, nOrd, msiz, nOsc, mMat(0:nOrd,nOsc), Graph2(m_max+1,m_max+1,nOsc)
integer(kind=iwp), intent(out) :: mInc(0:msiz,nOsc), mDec(0:msiz,nOsc)
integer(kind=iwp) :: i, iv, j
integer(kind=iwp), external :: iDetnr
integer(kind=iwp), allocatable :: iVec(:)

call mma_allocate(iVec,nOsc,label='iVec')

! Create mInc.

mInc(:,:) = -1
do i=0,msiz
  iVec(:) = mMat(i,:)
  do j=1,nOsc
    iVec(j) = iVec(j)+1
    mInc(i,j) = iDetnr(iVec,Graph2,nosc,m_max)
    iVec(j) = iVec(j)-1
  end do
end do

! Create mDec.

mDec(0,:) = -1
do i=1,msiz
  do j=1,nOsc
    if (mMat(i,j) /= 0) then
      iVec(:) = mMat(i,:)
      iVec(j) = iVec(j)-1
      mDec(i,j) = iDetnr(iVec,Graph2,nosc,m_max)
      do iv=1,nOsc
        iVec(iv) = iVec(j)+1
      end do
    else
      mDec(i,j) = -1
    end if
  end do
end do

call mma_deallocate(iVec)

return

end subroutine Mk_nIncDec
