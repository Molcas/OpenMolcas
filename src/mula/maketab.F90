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

implicit real*8(a-h,o-z)
integer mMat(0:maxord,nosc)
integer mInc(0:maxord,nosc), mDec(0:maxord,nosc)
logical equal
#include "WrkSpc.fh"

! Initialize.
maxOrd = 0
maxIncOrd = 0
do iv=0,maxord
  do jv=1,nosc
    mMat(iv,jv) = 0
    mInc(iv,jv) = 0
    mDec(iv,jv) = 0
  end do
end do
if (m_max == 0) return
if (nOsc == 1) then
  mTempDim = m_max
else
  mTempDim = (1-nOsc**(m_max+1))/(1-nOsc)
end if
lmTemp = mTempDim*nOsc
call GetMem('mTemp','Allo','INTE',ipmTemp,lmTemp)

call GetMem('row','Allo','INTE',iprow,nOsc)

!mTemp = 0
do iv=0,lmTemp-1
  iWork(ipmTemp+iv) = 0
end do
num = 1
istart_row = 0
irow = istart_row+1
call GetMem('unit','Allo','INTE',ipunit,nOsc*nOsc)
do i=1,nOsc
  do j=1,nOsc
    iWork(ipunit+i+nOsc*(j-1)-1) = 0
  end do
  iWork(ipunit+i+nOsc*(i-1)-1) = 1
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
        iWork(ipmTemp+irow+mTempDim*(jcol-1)) = mMat(istart_row,jcol)+iWork(ipunit+jrow+nOsc*(jcol-1)-1)
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
        if (iWork(ipmTemp+i+mTempDim*(l-1)) /= mMat(j,l)) then
          equal = .false.
        end if
      end do
      j = j+1
    end do
    if (.not. equal) then
      do jcol=1,nOsc
        mMat(mMat_row,jcol) = iWork(ipmTemp+i+mTempDim*(jcol-1))
      end do
      mMat_row = mMat_row+1
    end if
  end do
  num = mMat_row-istart_row
end do
maxOrd = mMat_row-1
call GetMem('mTemp','Free','INTE',ipmTemp,lmTemp)

! Create mInc.
maxIncOrd = maxOrd-num
do i=0,maxIncOrd
  do iv=1,nOsc
    iWork(iprow+iv-1) = mMat(i,iv)
  end do
  do j=1,nOsc
    iWork(iprow+j-1) = iWork(iprow+j-1)+1
    equal = .false.
    k = i+1
    do while ((.not. equal) .and. (k <= maxOrd))
      !row2 => mMat(k,:)
      equal = .true.
      do l=1,nOsc
        if (iWork(iprow+l-1) /= mMat(k,l)) then
          equal = .false.
        end if
      end do
      if (equal) then
        mInc(i,j) = k
      else
        k = k+1
      end if
    end do
    iWork(iprow+j-1) = iWork(iprow+j-1)-1
  end do
end do

! Create mDec.
do i=1,maxOrd
  do iv=1,nOsc
    iWork(iprow+iv-1) = mMat(i,iv)
  end do
  do j=1,nOsc
    if (iWork(iprow+j-1) > 0) then
      iWork(iprow+j-1) = iWork(iprow+j-1)-1
      equal = .false.
      k = 0
      do while ((.not. equal) .and. (k <= maxOrd))
        !row2 => mMat(k,:)
        equal = .true.
        do l=1,nOsc
          if (iWork(iprow+l-1) /= mMat(k,l)) then
            equal = .false.
          end if
        end do
        if (equal) then
          mDec(i,j) = k
        else
          k = k+1
        end if
      end do
      iWork(iprow+j-1) = iWork(iprow+j-1)+1
    else
      mDec(i,j) = 0
    end if
  end do
end do

call GetMem('row','Free','INTE',iprow,nOsc)
call GetMem('unit','Free','INTE',ipunit,nOsc*nOsc)

end subroutine MakeTab
