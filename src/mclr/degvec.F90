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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DEGVEC(VEC,NDIM,NDGVL,IDEG)
! A vector VEC is given with elements in ascending order
! group elements in degenerate pairs
!
!=======
! Input
!=======
! VEC : input vector
! NDIM : Number of elements in vec
!
!========
! Output
!========
! NDGVL : Number of degenerate values
! IDEG(I) : Number of elements in VEC with degenerate value I
!
! Jeppe Olsen, April 1990

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(in) :: VEC(*)
integer(kind=iwp), intent(in) :: NDIM
integer(kind=iwp), intent(out) :: NDGVL, IDEG(NDIM)
integer(kind=iwp) :: I, NDEG
real(kind=wp) :: XDGVL
real(kind=wp), parameter :: THRES = 1.0e-8_wp

!? write(u6,*) ' Input vector to DEGVEC'
!? call wrtmat(VEC,1,NDIM,1,NDIM)
XDGVL = VEC(1)
NDEG = 1
NDGVL = 0
do I=2,NDIM
  if (abs(VEC(I)-XDGVL) <= THRES) then
    NDEG = NDEG+1
  else
    NDGVL = NDGVL+1
    IDEG(NDGVL) = NDEG
    XDGVL = VEC(I)
    NDEG = 1
  end if
end do
! Last group
NDGVL = NDGVL+1
IDEG(NDGVL) = NDEG

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from DEGVEC'
write(u6,*) ' =================='
write(u6,*)
write(u6,*) ' Number of degenerate values ',NDGVL
write(u6,*) ' Degenerencies of each value'
call IWRTMA(IDEG,1,NDGVL,1,NDGVL)
#endif

return

end subroutine DEGVEC
