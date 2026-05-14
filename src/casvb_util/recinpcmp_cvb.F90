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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

function recinpcmp_cvb(ifield)

use casvb_global, only: recinp, recinp_old
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: recinpcmp_cvb
integer(kind=iwp), intent(in) :: ifield
integer(kind=iwp) :: i, ioff1, ioff2, joff1, joff2
logical(kind=iwp) :: done
real(kind=wp), allocatable :: tmp1(:), tmp2(:)
logical(kind=iwp), external :: valid_cvb ! ... Files/Hamiltonian available ...

if (.not. valid_cvb(recinp_old)) then
  recinpcmp_cvb = .true.
else
  call rdioff_cvb(ifield,recinp,ioff1)
  call rdioff_cvb(ifield+1,recinp,ioff2)
  call rdioff_cvb(ifield,recinp_old,joff1)
  call rdioff_cvb(ifield+1,recinp_old,joff2)
  if (ioff2-ioff1 /= joff2-joff1) then
    recinpcmp_cvb = .true.
  else
    call mma_allocate(tmp1,ioff2-ioff1,label='tmp1')
    call mma_allocate(tmp2,joff2-joff1,label='tmp2')
    call rdlow_cvb(tmp1,ioff2-ioff1,recinp,ioff1)
    call rdlow_cvb(tmp2,joff2-joff1,recinp_old,joff1)
    done = .false.
    do i=1,ioff2-ioff1
      if (tmp1(i) /= tmp2(i)) then
        recinpcmp_cvb = .true.
        done = .true.
        exit
      end if
    end do
    if (.not. done) recinpcmp_cvb = .false.
    call mma_deallocate(tmp1)
    call mma_deallocate(tmp2)
  end if
end if

return

end function recinpcmp_cvb
