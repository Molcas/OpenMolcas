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

subroutine dkh_geneu(n,m,xord,c,w,xl,xs,t1,t2,t3)
! Calculate the DKH unitary transformation truncated at [xord] order
!   U_{DKH}=U_{0}U_{1}U_{2}...U_{xord}
!   U_{k}=\sum_{i=0}^{xord/k}c_{i}W_{k}^{i}

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
! Input
!   n    dimension of matrix
!   m    =2*n
!   c    expansion coefficients of unitary transformation in terms of W
!   w    stored W matrices
! Output
!   xl   upper part
!   xs   lower part
integer(kind=iwp), intent(in) :: n, m, xord
real(kind=wp), intent(in) :: c(*), w(n,n,2,xord)
real(kind=wp), intent(out) :: xl(n,n), xs(n,n), t1(m,m), t2(m,m), t3(m,m)
integer(kind=iwp) :: i, j, k, iord

do iord=1,xord
  ! initial unit matrix
  call unitmat(t2,m)
  do k=1,xord/iord
    if (mod(k,2) == 1) then
      if (k == 1) then
        ! unitary transformation in right side, had {\dag} to original definition of U
        xs(:,:) = -w(:,:,1,iord)
      else
        call dmxma(n,'N','N',xl,w(1,1,1,iord),xs,-One)
      end if
      do i=1,n
        do j=1,n
          t2(j,i+n) = t2(j,i+n)+xs(j,i)*c(k)
          t2(j+n,i) = t2(j+n,i)-xs(i,j)*c(k)
        end do
      end do
    else
      call dmxma(n,'C','N',w(1,1,1,iord),xs,xl,One)
      do i=1,n
        do j=1,n
          t2(j+n,i+n) = t2(j+n,i+n)+xl(j,i)*c(k)
        end do
      end do
      call dmxma(n,'N','C',xs,w(1,1,1,iord),xl,One)
      do i=1,n
        do j=1,n
          t2(j,i) = t2(j,i)+xl(j,i)*c(k)
        end do
      end do
    end if
  end do
  if (iord == 1) then
    t1(:,:) = t2(:,:)
  else
    ! multiply U_{iord} in right side
    call dmxma(m,'N','N',t1,t2,t3,One)
    t1(:,:) = t3(:,:)
  end if
end do
do i=1,n
  do j=1,n
    xl(j,i) = t1(j,i)
    xs(j,i) = t1(j+n,i)
  end do
end do

return

end subroutine dkh_geneu
