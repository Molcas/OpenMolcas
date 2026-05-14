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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: R(3,3)
complex(kind=wp), intent(inout) :: S(3,n,n), M(3,n,n)
real(kind=wp), intent(out) :: mg(3,3)
logical(kind=iwp), intent(in) :: dbg
!real(kind=wp) :: g(3)
complex(kind=wp), allocatable :: Mt(:,:,:), St(:,:,:) !, Z(:,:)

!-----------------------------------------------------------------------
call mma_allocate(Mt,3,n,n,'Mt')
call mma_allocate(St,3,n,n,'St')
!call mma_allocate(Z,n,n,'Z')

! make a local backup of the data:
Mt(:,:,:) = M(:,:,:)
St(:,:,:) = S(:,:,:)
call unitmat(mg,3)

if (dbg) then
  call prMom('PA_prep_mom_exch, input S',St,n)
  call prMom('PA_prep_mom_exch, input M',Mt,n)
end if

! rotate the momentum using the R rotation matrix --
! to the local axes for a symmetric compound:
call rotmom2(St,n,R,S)
call rotmom2(Mt,n,R,M)
! back-up again:
!Mt(:,:,:) = M(:,:,:)
!St(:,:,:) = S(:,:,:)

!-----------------------------------------------------------------------
! experimental:
!! find local magnetic axes:
!call atens(M,n,g,mg,2)
!! rotate the momentum using the  mg  rotation matrix --
!! to the local magnetic axes:
!call rotmom2(St,n,mg,S)
!call rotmom2(Mt,n,mg,M)
!
!! find local pseudospin:
!call pseudospin(M,n,Z,3,1,1)
!if (dbg) call pa_prmat('PA_prep_mom_exch, Z:',Z,n)
!
!! Transform the moment into their local pseudospins
!call UTMU2(n,n,Z,S)
!call UTMU2(n,n,Z,M)
!if (dbg) then
!  call prMom('PA_prep_mom_exch, S:',S,n)
!  call prMom('PA_prep_mom_exch, M:',M,n)
!end if
!! back-up again:
!Mt(:,:,:) = M(:,:,:)
!St(:,:,:) = S(:,:,:)
!
!! rotate back the moment, so that we preserve the
!! original coordinate system of the computed molecule
!call rotmom(St,n,mg,S)
!call rotmom(Mt,n,mg,M)
!-----------------------------------------------------------------------

call mma_deallocate(Mt)
call mma_deallocate(St)
!call mma_deallocate(Z)

!-----------------------------------------------------------------------
! old preparation of the data for Lines exchange
!! rotate the moments to the general coordinate system
!
!SM_tmp(:,:,:) = SM(i1,:,1:n1,1:n1)
!r1(:,:) = rot(i1,j1,:,:)
!call rotmom2(SM_tmp,n1,r1,S1_tmp)
!SM_tmp(:,:,:) = SM(i2,:,1:n2,1:n2)
!r2(:,:) = rot(i2,j2,:,:)
!call rotmom2(SM_tmp,n2,r2,S2_tmp)
!call Lines_Exchange(Jex(lp),n1,n2,S1_tmp,S2_tmp,HTMP)
!HLIN1(lp,1:n1,1:n1,1:n2,1:n2) = HTMP(:,:,:,:)
!S1(:,1:n1,1:n1) = S1_tmp(:,:)
!S2(:,1:n2,1:n2) = S2_tmp(:,:)
!-----------------------------------------------------------------------

return

end subroutine prep_mom_exchange
