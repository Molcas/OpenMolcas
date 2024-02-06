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

subroutine transHam(n1,n2,rot1,rot2,MM1,MM2,typ1,typ2,H,HT,iopt)
! purpose:  transform the exchange Hamiltonian
!           matrices to the basis of their local pseudospins

use Constants, only: Zero, One, cZero, cOne
use Definitions, only: u6

implicit none
#include "stdalloc.fh"
! exchange basis of both sites
integer, intent(in) :: n1, n2, iopt
real(kind=8), intent(in) :: rot1(3,3)
real(kind=8), intent(in) :: rot2(3,3)
complex(kind=8), intent(in) :: MM1(3,n1,n1)
complex(kind=8), intent(in) :: MM2(3,n2,n2)
complex(kind=8), intent(in) :: H(n1,n1,n2,n2)
complex(kind=8), intent(out) :: HT(n1,n1,n2,n2)
character, intent(in) :: typ1, typ2
! local variables
integer :: i1, i2, j1, j2, i
real(kind=8), allocatable :: gt1(:), gt2(:)
real(kind=8), allocatable :: ax1(:,:), ax2(:,:)
complex(kind=8), allocatable :: M1(:,:,:), M2(:,:,:)
complex(kind=8), allocatable :: MR1(:,:,:), MR2(:,:,:)
complex(kind=8), allocatable :: Z1(:,:), Z2(:,:)
complex(kind=8), allocatable :: TMP1(:,:), TMP2(:,:)
complex(kind=8), allocatable :: HI(:,:,:,:) !HI(n1,n1,n2,n2)
logical :: DBG

DBG = .false.

!-----------------------------------------------------------------------
call mma_allocate(gt1,3,'gt1')
call mma_allocate(gt2,3,'gt2')
call mma_allocate(ax1,3,3,'ax1')
call mma_allocate(ax2,3,3,'ax2')
if (n1 > 0) then
  call mma_allocate(M1,3,n1,n1,'M1')
  call mma_allocate(MR1,3,n1,n1,'MR1')
  call mma_allocate(Z1,n1,n1,'Z1')
  call mma_allocate(TMP1,n1,n1,'TMP1')
end if
if (n2 > 0) then
  call mma_allocate(M2,3,n2,n2,'M2')
  call mma_allocate(MR2,3,n2,n2,'MR2')
  call mma_allocate(Z2,n2,n2,'Z2')
  call mma_allocate(TMP2,n2,n2,'TMP2')
end if
if ((n1 > 0) .and. (n2 > 0)) call mma_allocate(HI,n1,n1,n2,n2,'HI')
call zcopy_(n1*n1*n2*n2,[cZero],0,HT,1)
!-----------------------------------------------------------------------

! rotate the magnetic moments to their general coordinate system
call zcopy_(3*n1*n1,[cZero],0,M1,1)
call zcopy_(3*n2*n2,[cZero],0,M2,1)
call rotmom(MM1,n1,rot1,M1)
call rotmom(MM2,n2,rot2,M2)

if (iopt == 1) then ! local coordinate system
  call dcopy_(3,[Zero],0,gt1,1)
  call dcopy_(3,[Zero],0,gt2,1)
  call dcopy_(3*3,[Zero],0,ax1,1)
  call dcopy_(3*3,[Zero],0,ax2,1)
  if (DBG) then
    write(u6,'(A)') 'TRANSHAM:: local g tensors and axes:'
    call atens(M1,n1,gt1,ax1,2)
    call atens(M2,n2,gt2,ax2,2)
  else
    call atens(M1,n1,gt1,ax1,1)
    call atens(M2,n2,gt2,ax2,1)
  end if

else if (iopt == 2) then ! general coordinate system
  call dcopy_(3*3,[Zero],0,ax1,1)
  call dcopy_(3*3,[Zero],0,ax2,1)
  do i1=1,3
    ax1(i1,i1) = One
    ax2(i1,i1) = One
  end do
end if

!-----------------------------------------------------------------------
! rotate magnetic moments to their local magnetic axes:
call zcopy_(3*n1*n1,[cZero],0,MR1,1)
call zcopy_(3*n2*n2,[cZero],0,MR2,1)
if (((typ1 == 'A') .and. (typ2 == 'A')) .or. ((typ1 == 'B') .and. (typ2 == 'B')) .or. ((typ1 == 'B') .and. (typ2 == 'C')) .or. &
    ((typ1 == 'C') .and. (typ2 == 'B')) .or. ((typ1 == 'C') .and. (typ2 == 'C'))) then
  ! both sites have been computed ab initio
  call rotmom2(M1,n1,ax1,MR1)
  call rotmom2(M2,n2,ax2,MR2)
else if (((typ1 == 'A') .and. (typ2 == 'B')) .or. ((typ1 == 'A') .and. (typ2 == 'C'))) then
  ! site 2 is generated isotropic
  ! use axes of site 1
  call rotmom2(M1,n1,ax1,MR1)
  call rotmom2(M2,n2,ax1,MR2)
else if (((typ1 == 'B') .and. (typ2 == 'A')) .or. ((typ1 == 'C') .and. (typ2 == 'A'))) then
  ! site 1 is generated isotropic
  ! use axes of site 2
  call rotmom2(M1,n1,ax2,MR1)
  call rotmom2(M2,n2,ax2,MR2)
end if

if (DBG) then
  if (iopt == 1) then
    write(u6,'(A,i3)') 'site 1'
    call prMom('transHam:: magnetic moment, coordinate system of LOCAL main magnetic axes, site 1',MR1,n1)
    write(u6,'(A,i3)') 'site 2'
    call prMom('transHam:: magnetic moment, coordinate system of LOCAL main magnetic axes, site 2',MR2,n2)
  else
    write(u6,'(A,i3)') 'site 1'
    call prMom('transHam:: magnetic moment, coordinate system of GENERAL main magnetic axes, site 1',MR1,n1)
    write(u6,'(A,i3)') 'site 2'
    call prMom('transHam:: magnetic moment, coordinate system of GENERAL main magnetic axes, site 2',MR2,n2)
  end if
end if

!-----------------------------------------------------------------------
! find local pseudospin on each site:
call zcopy_(n1*n1,[cZero],0,Z1,1)
call zcopy_(n2*n2,[cZero],0,Z2,1)
call pseudospin(MR1,n1,Z1,3,1,1)
call pseudospin(MR2,n2,Z2,3,1,1)

if (DBG) then
  call pa_prmat('Matrix Z1',Z1,n1)
  call zgemm_('C','N',n1,n1,n1,cOne,Z1,n1,Z1,n1,cZero,TMP1,n1)
  write(u6,'(A)') 'transHam:  Verify unitarity of Z1'
  do i=1,n1
    write(u6,'(A,i2,A,i2,A,2ES20.10)') 'conjg(Z1)*Z1:  (',i,',',i,')=',TMP1(i,i)
  end do
  call pa_prmat('Matrix Z2',Z2,n2)
  call zgemm_('C','N',n2,n2,n2,cOne,Z2,n2,Z2,n2,cZero,TMP2,n2)
  write(u6,'(A)') 'transHam:  Verify unitarity of Z2'
  do i=1,n2
    write(u6,'(A,i2,A,i2,A,2ES20.10)') 'conjg(Z2)*Z2:  (',i,',',i,')=',TMP2(i,i)
  end do
  call UTMU2(n1,n1,Z1,MR1)
  call UTMU2(n2,n2,Z2,MR2)
  call prMom('transHam:: magnetic moment, coordinate MR1, site 1',MR1,n1)
  call prMom('transHam:: magnetic moment, coordinate MR2, site 2',MR2,n2)
end if

! save a local copy:
call zcopy_(n1*n1*n2*n2,H,1,HI,1)

do i1=1,n1
  do j1=1,n1
    do i2=1,n2
      do j2=1,n2
        !if (DBG .and. (abs(H(i1,j1,i2,j2)) > 0.5e-13_wp)) then
        write(u6,'(A,4(i2,A),2ES22.14)') 'H (',i1,',',j1,',',i2,',',j2,')=',H(i1,j1,i2,j2)
        !end if
      end do
    end do
  end do
end do

! transform the Hamiltonian to local pseudospins:
do i2=1,n2
  do j2=1,n2
    call zgemm_('C','N',n1,n1,n1,cOne,Z1,n1,HI(:,:,i2,j2),n1,cZero,TMP1,n1)
    call zgemm_('N','N',n1,n1,n1,cOne,TMP1,n1,Z1,n1,cZero,HI(:,:,i2,j2),n1)
  end do
end do

do i1=1,n1
  do j1=1,n1
    call zgemm_('C','N',n2,n2,n2,cOne,Z2,n2,HI(i1,j1,:,:),n2,cZero,TMP2,n2)
    call zgemm_('N','N',n2,n2,n2,cOne,TMP2,n2,Z2,n2,cZero,HI(i1,j1,:,:),n2)
  end do
end do

call zcopy_(n1*n1*n2*n2,[cZero],0,HT,1)
call zcopy_(n1*n1*n2*n2,HI,1,HT,1)

do i1=1,n1
  do j1=1,n1
    do i2=1,n2
      do j2=1,n2
        !if (DBG .and. (abs(HT(i1,j1,i2,j2)) > 0.5e-13_wp)) then
        write(u6,'(A,4(i2,A),2ES22.14)') 'HT(',i1,',',j1,',',i2,',',j2,')=',HT(i1,j1,i2,j2)
        !end if
      end do
    end do
  end do
end do

!-----------------------------------------------------------------------
call mma_deallocate(gt1)
call mma_deallocate(gt2)
call mma_deallocate(ax1)
call mma_deallocate(ax2)
if (n1 > 0) then
  call mma_deallocate(M1)
  call mma_deallocate(MR1)
  call mma_deallocate(Z1)
  call mma_deallocate(TMP1)
end if
if (n2 > 0) then
  call mma_deallocate(M2)
  call mma_deallocate(MR2)
  call mma_deallocate(Z2)
  call mma_deallocate(TMP2)
end if
if ((n1 > 0) .and. (n2 > 0)) call mma_deallocate(HI)

return

end subroutine transHam
