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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine MakeTab(m_max,maxOrd,maxIncOrd,mMat,mInc,mDec,nOsc)
!  Purpose:
!    Create tables used in FCval.
!
!  Input:
!    nOsc      : Integer - the the number of oscillators.
!    m_max     : Integer - the maximum sum of the quantum numbers.
!    maxOrd    : Integer - number of rows in mMat.
!    nTabDim   : Integer
!    Osc_Shift : Integer array
!
!  Output:
!    mMat      : Two dimensional integer array
!    mInc      : Two dimensional integer array
!    mDec      : Two dimensional integer array
!
!  Calls:
!    none
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate

implicit real*8(a-h,o-z)
integer mMat(0:maxord,nosc)
integer mInc(0:maxord,nosc), mDec(0:maxord,nosc)
logical equal
integer, allocatable :: mTemp(:,:), row(:), unit(:,:)

! Initialize.
maxOrd = 0
maxIncOrd = 0
mMat(:,:) = 0
mInc(:,:) = 0
mDec(:,:) = 0
if (m_max == 0) return
if (nOsc == 1) then
  mTempDim = m_max
else
  mTempDim = (1-nOsc**(m_max+1))/(1-nOsc)
end if
call mma_allocate(mTemp,[0,mTempDim],[1,nOsc],label='mTemp')

call mma_allocate(row,nOsc,label='row')

mTemp(:,:) = 0
num = 1
istart_row = 0
irow = istart_row+1
call mma_allocate(unit,nOsc,nOsc,label='unit')
unit(:,:) = 0
do i=1,nOsc
  unit(i,i) = 1
end do

! Create table mMat.
mMat_row = 1
do m=1,m_max
  irow = 1
  numtemp = 0
  ! Produce all combinations for a given total sum = m_max.
  do n=1,num
    do jrow=1,nOsc
      do jcol=1,nOsc
        mTemp(irow,jcol) = mMat(istart_row,jcol)+unit(jrow,jcol)
      end do
      irow = irow+1
      numtemp = numtemp+1
    end do
    istart_row = istart_row+1
  end do
  num = numtemp
  ! Remove all entries which occur more than once.
  do i=1,num
    !row1 => mTemp(i,:)
    j = istart_row
    jmax = mMat_row
    equal = .false.
    do while ((j < jmax) .and. (.not. equal))
      !row2 => mMat(j,:)
      equal = .true.
      do l=1,nOsc
        !if (row1(l) /= row2(l)) then
        if (mTemp(i,l) /= mMat(j,l)) then
          equal = .false.
        end if
      end do
      j = j+1
    end do
    if (.not. equal) then
      do jcol=1,nOsc
        mMat(mMat_row,jcol) = mTemp(i,jcol)
      end do
      mMat_row = mMat_row+1
    end if
  end do
  num = mMat_row-istart_row
end do
maxOrd = mMat_row-1
call mma_deallocate(mTemp)

! Create mInc.
maxIncOrd = maxOrd-num
do i=0,maxIncOrd
  row(:) = mMat(i,:)
  do j=1,nOsc
    row(j) = row(j)+1
    equal = .false.
    k = i+1
    do while ((.not. equal) .and. (k <= maxOrd))
      !row2 => mMat(k,:)
      equal = .true.
      do l=1,nOsc
        if (row(l) /= mMat(k,l)) then
          equal = .false.
        end if
      end do
      if (equal) then
        mInc(i,j) = k
      else
        k = k+1
      end if
    end do
    row(j) = row(j)-1
  end do
end do

! Create mDec.
do i=1,maxOrd
  row(:) = mMat(i,:)
  do j=1,nOsc
    if (row(j) > 0) then
      row(j) = row(j)-1
      equal = .false.
      k = 0
      do while ((.not. equal) .and. (k <= maxOrd))
        !row2 => mMat(k,:)
        equal = .true.
        do l=1,nOsc
          if (row(l) /= mMat(k,l)) then
            equal = .false.
          end if
        end do
        if (equal) then
          mDec(i,j) = k
        else
          k = k+1
        end if
      end do
      row(j) = row(j)+1
    else
      mDec(i,j) = 0
    end if
  end do
end do

call mma_deallocate(row)
call mma_deallocate(unit)

end subroutine MakeTab
