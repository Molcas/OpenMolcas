!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine MergeBS(z1,n1,z2,n2,z,n,RatioThres,iDominantSet)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n1, n2, iDominantSet
real(kind=wp), intent(in) :: z1(n1), z2(n2), RatioThres
real(kind=wp), intent(_OUT_) :: z(*)
integer(kind=iwp), intent(out) :: n
integer(kind=iwp) :: i, i1, i2, idum, iSetPrev, iSetThis, j, mPrim
real(kind=wp) :: ratio
integer(kind=iwp), allocatable :: ix1(:), ix2(:)
#ifdef _DEBUGPRINT_
#define _TEST_ .true.
#else
#define _TEST_ .false.
#endif
logical(kind=iwp), parameter :: IfTest = _TEST_

mPrim = n1+n2
call mma_allocate(ix1,mPrim,label='ix1')
call mma_allocate(ix2,mPrim,label='ix2')

iSetPrev = 0

do i=1,mPrim
  ix1(i) = i
  ix2(i) = i
end do

do i=1,n1-1
  do j=i+1,n1
    if (z1(ix1(i)) < z1(ix1(j))) then
      idum = ix1(i)
      ix1(i) = ix1(j)
      ix1(j) = idum
    end if
  end do
end do

do i=1,n2-1
  do j=i+1,n2
    if (z2(ix2(i)) < z2(ix2(j))) then
      idum = ix2(i)
      ix2(i) = ix2(j)
      ix2(j) = idum
    end if
  end do
end do

if (IfTest) then
  write(u6,'(A)')
  write(u6,'(4f20.4)') (z1(ix1(i)),i=1,n1)
  write(u6,'(A)')
  write(u6,'(4f20.4)') (z2(ix2(i)),i=1,n2)
end if

i = 0
i1 = 1
i2 = 1
do while ((i1 <= n1) .or. (i2 <= n2))
  i = i+1
  if (i > mPrim) then
    call WarningMessage(2,'Error in MergeBS')
    write(u6,*) ' MergeBS: i > mPrim',i,mPrim
    write(u6,*) ' MergeBS: raise mPrim and recompile'
    call Abend()
  end if
  if (i1 > n1) then
    z(i) = z2(ix2(i2))
    iSetThis = 2
    i2 = i2+1
  else if (i2 > n2) then
    z(i) = z1(ix1(i1))
    iSetThis = 1
    i1 = i1+1
  else if (z1(ix1(i1)) > z2(ix2(i2))) then
    z(i) = z1(ix1(i1))
    iSetThis = 1
    i1 = i1+1
  else
    z(i) = z2(ix2(i2))
    iSetThis = 2
    i2 = i2+1
  end if

  if (i == 1) then
    iSetPrev = iSetThis
  else
    ratio = z(i-1)/z(i)
    if (ratio >= RatioThres) then
      iSetPrev = iSetThis
    else
      if (iSetThis /= iDominantSet) then
        i = i-1
      else if (iSetThis == iDominantSet) then
        if (iSetPrev /= iDominantSet) then
          z(i-1) = z(i)
          i = i-1
        end if
        iSetPrev = iSetThis
      end if
    end if
  end if

end do

call mma_deallocate(ix1)
call mma_deallocate(ix2)

n = i
if (IfTest) then
  write(u6,'(I4)') n
  write(u6,'(4f20.4)') (z(i),i=1,n)
end if

end subroutine MergeBS
