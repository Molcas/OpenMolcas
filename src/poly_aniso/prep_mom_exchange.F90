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

subroutine prep_mom_exchange(n,R,S,M,mg,dbg)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: n
real(kind=8), intent(in) :: R(3,3)
real(kind=8), intent(out) :: mg(3,3)
complex(kind=8), intent(inout) :: S(3,n,n), M(3,n,n)
logical :: dbg
! local data:
integer :: i
complex(kind=8), allocatable :: Mt(:,:,:), St(:,:,:)
!real(kind=8) :: g(3)
!complex(kind=8), allocatable :: Z(:,:)

!-----------------------------------------------------------------------
call mma_allocate(Mt,3,n,n,'Mt')
call mma_allocate(St,3,n,n,'St')
!call mma_allocate(Z,n,n,'Z')
!call dcopy_(3,[0.0_wp],0,g,1)
call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,Mt,1)
call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,St,1)
call dcopy_(3*3,[0.0_wp],0,mg,1)

! make a local backup of the data:
call zcopy_(3*n*n,M,1,Mt,1)
call zcopy_(3*n*n,S,1,St,1)
do i=1,3
  mg(i,i) = 1.0_wp
end do

if (dbg) call prMom('PA_prep_mom_exch, input S',St,n)
if (dbg) call prMom('PA_prep_mom_exch, input M',Mt,n)

! rotate the momentum using the R rotation matrix --
! to the local axes for a symmetric compound:
call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
call rotmom2(St,n,R,S)
call rotmom2(Mt,n,R,M)
! back-up again:
!call zcopy_(3*n*n,M,1,Mt,1)
!call zcopy_(3*n*n,S,1,St,1)

!-----------------------------------------------------------------------
! experimental:
!if (.false.) then
!  ! find local magnetic axes:
!  call atens(M,n,g,mg,2)
!  ! rotate the momentum using the  mg  rotation matrix --
!  ! to the local magnetic axes:
!  call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
!  call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
!  call rotmom2(St,n,mg,S)
!  call rotmom2(Mt,n,mg,M)
!
!  ! find local pseudospin:
!  call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Z,1)
!  call pseudospin(M,n,Z,3,1,1)
!  if (dbg) call pa_prmat('PA_prep_mom_exch, Z:',Z,n)
!
!  ! Transform the moment into their local pseudospins
!  call UTMU2(n,n,Z,S)
!  call UTMU2(n,n,Z,M)
!  if (dbg) call prMom('PA_prep_mom_exch, S:',S,n)
!  if (dbg) call prMom('PA_prep_mom_exch, M:',M,n)
!  ! back-up again:
!  call zcopy_(3*n*n,M,1,Mt,1)
!  call zcopy_(3*n*n,S,1,St,1)
!
!  ! rotate back the moment, so that we preserve the
!  ! original coordinate system of the computed molecule
!  call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
!  call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
!  call rotmom(St,n,mg,S)
!  call rotmom(Mt,n,mg,M)
!end if
!-----------------------------------------------------------------------

call mma_deallocate(Mt)
call mma_deallocate(St)
!call mma_deallocate(Z)

!-----------------------------------------------------------------------
! old preparation of the data for Lines exchange
!! rotate the moments to the general coordinate system
!call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S1,1)
!call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S2,1)
!
!call rotmom2(SM(i1,1:3,1:n1,1:n1),n1,rot(i1,j1,1:3,1:3),S1(1:3,1:n1,1:n1))
!call rotmom2(SM(i2,1:3,1:n2,1:n2),n2,rot(i2,j2,1:3,1:3),S2(1:3,1:n2,1:n2))
!call Lines_Exchange(Jex(lp),n1,n2,S1(1:3,1:n1,1:n1),S2(1:3,1:n2,1:n2),HLIN1(lp,1:n1,1:n1,1:n2,1:n2))
!-----------------------------------------------------------------------

return

end subroutine prep_mom_exchange
