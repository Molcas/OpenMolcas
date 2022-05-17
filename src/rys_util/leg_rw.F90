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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

! Compute and store roots and weights for a shifted Legendre quadrature,
! used for boot-strapping Rys roots and weights.
! Different sets of roots and weights are computed

module Leg_RW

implicit none
integer, dimension(11), parameter :: naux = [30,35,40,45,50,55,60,65,70,75,300]
real*8, dimension(:,:), allocatable :: Leg_r, Leg_w

contains

subroutine SetAux(eps)

  real*8, intent(In) :: eps
  integer, parameter :: nquad = size(naux)
  real*8, dimension(:), allocatable :: a, b
  integer :: maux, i, j, Err
# include "stdalloc.fh"
# include "real.fh"

  if (allocated(Leg_r)) return
  maux = maxval(naux)
  call mma_allocate(Leg_r,maux,nquad,label='Leg_r')
  call mma_allocate(Leg_w,maux,nquad,label='Leg_w')
  call mma_allocate(a,maux)
  call mma_allocate(b,maux)
  do j=1,nquad
    do i=1,naux(j)
      a(i) = Half
      if (i == 1) then
        b(1) = One
      else
        b(i) = Quart/(Four-One/(i-1)**2)
      end if
    end do
    call GaussQuad(naux(j),a,b,eps,Leg_r(1,j),Leg_w(1,j),Err)
    if (Err /= 0) then
      write(6,*) Err
      call WarningMessage(2,'Error in GaussQuad')
      call AbEnd()
    end if
    do i=1,naux(j)
      Leg_r(i,j) = Leg_r(i,j)*Leg_r(i,j)
    end do
  end do
  call mma_deallocate(a)
  call mma_deallocate(b)

end subroutine SetAux

subroutine UnSetAux()

# include "stdalloc.fh"

  if (allocated(Leg_r)) call mma_deallocate(Leg_r)
  if (allocated(Leg_w)) call mma_deallocate(Leg_w)

end subroutine UnSetAux

end module Leg_RW
