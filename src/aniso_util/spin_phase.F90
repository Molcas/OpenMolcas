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

subroutine SPIN_PHASE(MM,dim,Zinp,Zout)
! The RASSI program gives a random phase to the spin-orbit functions.
!
! This routine performs a simple check with the obtained spin functions,
! in order to determine the phase of the spin functions.
! If the phase is not the same, Then the spin functions will be multiplied
! with the correspondind coefficient that sets the same phase to all spin
! eigenfunctions

use Constants, only: Zero, One, cZero, cOne
use Definitions, only: wp, u6

implicit none
#include "stdalloc.fh"
integer, intent(in) :: dim
complex(kind=8), intent(in) :: mm(3,dim,dim)
complex(kind=8), intent(in) :: Zinp(dim,dim)
complex(kind=8), intent(out) :: Zout(dim,dim)
!-----------------------------------------------------------------------
integer :: i, j, i1, i2, l
real(kind=8), allocatable :: rxr(:) !dim)
real(kind=8), allocatable :: rxi(:) !dim)
complex(kind=8), allocatable :: r(:) !(dim)
complex(kind=8), allocatable :: phs(:,:,:)  !3,dim,dim)
complex(kind=8), allocatable :: tmp(:,:) !dim,dim
logical :: dbg

dbg = .false.

call mma_allocate(rxr,dim,'rxr')
call mma_allocate(rxi,dim,'rxi')
call mma_allocate(r,dim,'r')
call mma_allocate(phs,3,dim,dim,'phs')
call mma_allocate(tmp,dim,dim,'tmp')
!-----------------------------------------------------------------------
phs(:,:,:) = cZero
r(:) = czero
rxr(:) = Zero
rxi(:) = Zero
rxi(1) = Zero
r(1) = cOne

do i=1,dim-1
  j = i+1
  Zout(:,1) = Zinp(:,1)

  do i1=1,dim
    do i2=1,dim
      phs(1,i,j) = phs(1,i,j)+MM(1,i1,i2)*conjg(Zout(i1,i))*Zinp(i2,j)
    end do
  end do

  if (abs(phs(1,i,j)) > 1.0e-14_wp) then
    rxr(j) = real(phs(1,i,j))/abs(phs(1,i,j))
    rxi(j) = aimag(phs(1,i,j))/abs(phs(1,i,j))
  else
    rxr(j) = One
    rxi(j) = Zero
  end if

  r(j) = cmplx(rxr(j),rxi(j),kind=wp)

  Zout(:,j) = conjg(r(j))*Zinp(:,j)

  if (dbg) write(u6,'(A,i2,A,2ES24.14)') 'SPIN-PHASE: R(',j,') = ',conjg(r(j))
end do ! i

call zgemm_('C','N',dim,dim,dim,cOne,Zout,dim,mm(1,:,:),dim,cZero,TMP,dim)
call zgemm_('N','N',dim,dim,dim,cOne,TMP,dim,Zout,dim,cZero,phs(1,:,:),dim)
! convention:
!    mX(i,i+1) => Real, negative
!    mY(i,i+1) => imag, positive
!    mZ(i,i)   => diagonal
do i=1,dim-1,2
  j = i+1
  if (real(phs(1,i,j)) > Zero) Zout(:,j) = -Zout(:,j)
end do

if (dbg) then

  do l=1,3
    call zgemm_('C','N',dim,dim,dim,cOne,Zout,dim,mm(l,:,:),dim,cZero,TMP,dim)
    call zgemm_('N','N',dim,dim,dim,cOne,TMP,dim,Zout,dim,cZero,phs(l,:,:),dim)
  end do

  do i=1,dim
    do j=1,dim
      write(u6,'(a,i2,a,i2,a,2ES24.14)') 'SPIN-PHASE:  Zout(',i,',',j,') = ',Zout(i,j)
    end do
  end do

  write(u6,'(//)')
  do i=1,dim
    do j=1,dim
      if ((j == i-1) .or. (j == i+1)) write(u6,'(A,i2,A,i2,A, 3(2ES24.14,3x))') 'SPIN-PHASE: PHS(',i,',',j,') = (x,y,z) =', &
                                                                                (phs(l,i,j),l=1,3)
    end do
  end do
end if

!-----------------------------------------------------------------------
call mma_deallocate(rxr)
call mma_deallocate(rxi)
call mma_deallocate(r)
call mma_deallocate(phs)
call mma_deallocate(tmp)

return

end subroutine SPIN_PHASE
