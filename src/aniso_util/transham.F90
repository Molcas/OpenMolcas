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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n1, n2, iopt
real(kind=wp), intent(in) :: rot1(3,3), rot2(3,3)
complex(kind=wp), intent(in) :: MM1(3,n1,n1), MM2(3,n2,n2), H(n1,n1,n2,n2)
character, intent(in) :: typ1, typ2
complex(kind=wp), intent(out) :: HT(n1,n1,n2,n2)
integer(kind=iwp) :: i1, i2, iprint, j1, j2
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp), allocatable :: ax1(:,:), ax2(:,:), gt1(:), gt2(:)
complex(kind=wp), allocatable :: HI(:,:,:,:), HTMP(:,:), M1(:,:,:), M2(:,:,:), MR1(:,:,:), MR2(:,:,:), TMP1(:,:), TMP2(:,:), &
                                 Z1(:,:), Z2(:,:)

!-----------------------------------------------------------------------
call mma_allocate(gt1,3,'gt1')
call mma_allocate(gt2,3,'gt2')
call mma_allocate(ax1,3,3,'ax1')
call mma_allocate(ax2,3,3,'ax2')
call mma_allocate(M1,3,n1,n1,'M1')
call mma_allocate(MR1,3,n1,n1,'MR1')
call mma_allocate(Z1,n1,n1,'Z1')
call mma_allocate(TMP1,n1,n1,'TMP1')
call mma_allocate(M2,3,n2,n2,'M2')
call mma_allocate(MR2,3,n2,n2,'MR2')
call mma_allocate(Z2,n2,n2,'Z2')
call mma_allocate(TMP2,n2,n2,'TMP2')
call mma_allocate(HI,n1,n1,n2,n2,'HI')
HT(:,:,:,:) = cZero
!-----------------------------------------------------------------------

! rotate the magnetic moments to their general coordinate system

call rotmom(MM1,n1,rot1,M1)
call rotmom(MM2,n2,rot2,M2)

if (iopt == 1) then ! local coordinate system
  iprint = 1
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'TRANSHAM:: local g tensors and axes:'
  iprint = 2
# endif
  call atens(M1,n1,gt1,ax1,iprint)
  call atens(M2,n2,gt2,ax2,iprint)

else if (iopt == 2) then ! general coordinate system
  call unitmat(ax1,3)
  call unitmat(ax2,3)
end if

!-----------------------------------------------------------------------
! rotate magnetic moments to their local magnetic axes:
MR1(:,:,:) = cZero
MR2(:,:,:) = cZero
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

#ifdef _DEBUGPRINT_
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
#endif

!-----------------------------------------------------------------------
! find local pseudospin on each site:
call pseudospin(MR1,n1,Z1,3,1,1)
call pseudospin(MR2,n2,Z2,3,1,1)

#ifdef _DEBUGPRINT_
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
#endif

! save a local copy:
HI(:,:,:,:) = H(:,:,:,:)

do i1=1,n1
  do j1=1,n1
    do i2=1,n2
      do j2=1,n2
#       ifdef _DEBUGPRINT_
        if (abs(H(i1,j1,i2,j2)) > 0.5e-13_wp) write(u6,'(A,4(i2,A),2ES22.14)') 'H (',i1,',',j1,',',i2,',',j2,')=',H(i1,j1,i2,j2)
#       else
        write(u6,'(A,4(i2,A),2ES22.14)') 'H (',i1,',',j1,',',i2,',',j2,')=',H(i1,j1,i2,j2)
#       endif
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

call mma_allocate(HTMP,n2,n2,label='HTMP')
do i1=1,n1
  do j1=1,n1
    HTMP(:,:) = HI(i1,j1,:,:)
    call zgemm_('C','N',n2,n2,n2,cOne,Z2,n2,HTMP,n2,cZero,TMP2,n2)
    call zgemm_('N','N',n2,n2,n2,cOne,TMP2,n2,Z2,n2,cZero,HTMP,n2)
    HI(i1,j1,:,:) = HTMP(:,:)
  end do
end do
call mma_deallocate(HTMP)

HT(:,:,:,:) = HI(:,:,:,:)

do i1=1,n1
  do j1=1,n1
    do i2=1,n2
      do j2=1,n2
#       ifdef _DEBUGPRINT_
        if (abs(HT(i1,j1,i2,j2)) > 0.5e-13_wp) write(u6,'(A,4(i2,A),2ES22.14)') 'HT(',i1,',',j1,',',i2,',',j2,')=',HT(i1,j1,i2,j2)
#       else
        write(u6,'(A,4(i2,A),2ES22.14)') 'HT(',i1,',',j1,',',i2,',',j2,')=',HT(i1,j1,i2,j2)
#       endif
      end do
    end do
  end do
end do

!-----------------------------------------------------------------------
call mma_deallocate(gt1)
call mma_deallocate(gt2)
call mma_deallocate(ax1)
call mma_deallocate(ax2)
call mma_deallocate(M1)
call mma_deallocate(MR1)
call mma_deallocate(Z1)
call mma_deallocate(TMP1)
call mma_deallocate(M2)
call mma_deallocate(MR2)
call mma_deallocate(Z2)
call mma_deallocate(TMP2)
call mma_deallocate(HI)

return

end subroutine transHam
