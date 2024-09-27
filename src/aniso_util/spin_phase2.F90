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

subroutine SPIN_PHASE2(MM,d,Zinp,Zout)
! The RASSI program gives a random phase to the spin-orbit functions.
!
! This routine performs a simple check with the obtained spin functions,
! in order to determine the phase of the spin functions.
! If the phase is not the same, Then the spin functions will be multiplied
! with the correspondind coefficient that sets the same phase to all spin
! eigenfunctions

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, cZero, cOne
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: d
complex(kind=wp), intent(in) :: mm(3,d,d), Zinp(d,d)
complex(kind=wp), intent(out) :: Zout(d,d)
integer(kind=iwp) :: i, j
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: l
#endif
complex(kind=wp) :: t
real(kind=wp), allocatable :: rxi(:), rxr(:)
complex(kind=wp), allocatable :: mm_tmp(:,:), phs(:,:,:), r(:), tmp(:,:)

call mma_allocate(rxr,d,'rxr')
call mma_allocate(rxi,d,'rxi')
call mma_allocate(r,d,'r')
call mma_allocate(phs,d,d,3,'phs')
call mma_allocate(tmp,d,d,'tmp')
call mma_allocate(mm_tmp,d,d,'mm_tmp')
!-----------------------------------------------------------------------
r(:) = cZero
rxr(:) = Zero
rxi(:) = Zero
rxr(1) = One
r(1) = cOne

! compute magnetic moment, X
mm_tmp(:,:) = mm(1,:,:)
call zgemm_('C','N',d,d,d,cOne,Zinp,d,mm_tmp,d,cZero,TMP,d)
call zgemm_('N','N',d,d,d,cOne,TMP,d,Zinp,d,cZero,phs(:,:,1),d)

do i=1,d-1
  j = i+1

  if (abs(phs(i,j,1)) > 1.0e-14_wp) then
    rxr(j) = real(phs(i,j,1))/abs(phs(i,j,1))
    rxi(j) = aimag(phs(i,j,1))/abs(phs(i,j,1))
  else
    rxr(j) = One
    rxi(j) = Zero
  end if

  r(j) = cmplx(rxr(j),-rxi(j),kind=wp)

# ifdef _DEBUGPRINT_
  write(u6,'(A,i2,A,2ES24.14)') 'SPIN-PHASE: R(',j,') = ',r(j)*r(i)
# endif
end do

t = cOne
do j=1,d
  t = t*r(j)
  Zout(:,j) = t*Zinp(:,j)
end do

! compute the momentum using the ZOUT functions:
call zgemm_('C','N',d,d,d,cOne,Zout,d,mm_tmp,d,cZero,TMP,d)
call zgemm_('N','N',d,d,d,cOne,TMP,d,Zout,d,cZero,phs(:,:,1),d)
! convention:
!    mX(i,i+1) => Real, negative
!    mY(i,i+1) => imag, positive
!    mZ(i,i)   => diagonal
do i=1,d-1,2
  j = i+1
  if (real(phs(i,j,1)) > Zero) Zout(:,j) = -Zout(:,j)
end do

#ifdef _DEBUGPRINT_
do l=1,3
  mm_tmp(:,:) = mm(l,:,:)
  call ZGEMM_('C','N',d,d,d,cOne,Zout,d,mm_tmp,d,cZero,TMP,d)
  call ZGEMM_('N','N',d,d,d,cOne,TMP,d,Zout,d,cZero,phs(:,:,l),d)
end do

do i=1,d
  do j=1,d
    write(u6,'(a,i2,a,i2,a,2ES24.14)') 'SPIN-PHASE:  Zout(',i,',',j,') = ',Zout(i,j)
  end do
end do

write(u6,'(//)')
do i=1,d
  do j=1,d
    if ((j == i-1) .or. (j == i+1)) &
      write(u6,'(A,i2,A,i2,A, 3(2ES24.14,3x))') 'SPIN-PHASE: PHS(',i,',',j,') = (x,y,z) =',(phs(i,j,l),l=1,3)
  end do
end do
#endif

!-----------------------------------------------------------------------
call mma_deallocate(rxr)
call mma_deallocate(rxi)
call mma_deallocate(r)
call mma_deallocate(phs)
call mma_deallocate(tmp)
call mma_deallocate(mm_tmp)

return

end subroutine SPIN_PHASE2
