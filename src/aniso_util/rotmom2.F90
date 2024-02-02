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

subroutine rotmom2(MOM,N,R,MOMR)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N
real(kind=8), intent(in) :: R(3,3) !rotation matrix
complex(kind=8), intent(in) :: MOM(3,N,N) !initial momentum matrix
! rotated momentum matrix
complex(kind=8), intent(out) :: MOMR(3,N,N)
! local variables
integer i, j, l, k
complex(kind=8) :: RC(3,3)

! rotate the matrix
call zcopy_(3*N*N,[(0.0_wp,0.0_wp)],0,MOMR,1)

do l=1,3
  do k=1,3
    RC(l,k) = (0.0_wp,0.0_wp)
    RC(l,k) = cmplx(R(l,k),0.0d0,wp)
  end do
end do

do i=1,N
  do j=1,N
    do l=1,3
      do k=1,3
        MOMR(l,i,j) = MOMR(l,i,j)+RC(k,l)*MOM(k,i,j)
      end do
    end do
  end do
end do

return

end subroutine rotmom2
