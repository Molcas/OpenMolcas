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

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
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
call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,phs,1)
call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,tmp,1)
call zcopy_(dim,[(0.0_wp,0.0_wp)],0,r,1)
call dcopy_(dim,[0.0_wp],0,rxr,1)
call dcopy_(dim,[0.0_wp],0,rxi,1)
rxr(1) = 1.0_wp
rxi(1) = 0.0_wp

do i=1,dim-1
  j = i+1
  r(j) = (0.0_wp,0.0_wp)
  do i1=1,dim
    Zout(i1,1) = Zinp(i1,1)
  end do

  do i1=1,dim
    do i2=1,dim
      phs(1,i,j) = phs(1,i,j)+MM(1,i1,i2)*conjg(Zout(i1,i))*Zinp(i2,j)
    end do
  end do

  if (abs(phs(1,i,j)) > 1.0e-14_wp) then
    rxr(j) = dble(phs(1,i,j))/abs(phs(1,i,j))
    rxi(j) = aimag(phs(1,i,j))/abs(phs(1,i,j))
  else
    rxr(j) = 1.0_wp
    rxi(j) = 0.0_wp
  end if
  r(1) = (1.0_wp,0.0_wp)

  ! kind=8, complex double precision
  r(j) = cmplx(rxr(j),rxi(j),kind=8)

  do i1=1,dim
    Zout(i1,j) = conjg(r(j))*Zinp(i1,j)
  end do

  if (dbg) write(6,'(A,i2,A,2ES24.14)') 'SPIN-PHASE: R(',j,') = ',conjg(r(j))
end do ! i

call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,phs,1)
call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,tmp,1)
call zgemm_('C','N',dim,dim,dim,(1.0_wp,0.0_wp),Zout(1:dim,1:dim),dim,mm(1,1:dim,1:dim),dim,(0.0_wp,0.0_wp),TMP(1:dim,1:dim),dim)
call zgemm_('N','N',dim,dim,dim,(1.0_wp,0.0_wp),TMP(1:dim,1:dim),dim,Zout(1:dim,1:dim),dim,(0.0_wp,0.0_wp),phs(1,1:dim,1:dim),dim)
! convention:
!    mX(i,i+1) => Real, negative
!    mY(i,i+1) => imag, positive
!    mZ(i,i)   => diagonal
do i=1,dim-1,2
  j = i+1
  if (dble(phs(1,i,j)) > 0.0_wp) then
    do i1=1,dim
      Zout(i1,j) = -Zout(i1,j)
    end do
  end if
end do

if (dbg) then

  call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,phs,1)
  do l=1,3
    call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,tmp,1)
    call zgemm_('C','N',dim,dim,dim,(1.0_wp,0.0_wp),Zout(1:dim,1:dim),dim,mm(l,1:dim,1:dim),dim,(0.0_wp,0.0_wp),TMP(1:dim,1:dim), &
                dim)
    call zgemm_('N','N',dim,dim,dim,(1.0_wp,0.0_wp),TMP(1:dim,1:dim),dim,Zout(1:dim,1:dim),dim,(0.0_wp,0.0_wp),phs(l,1:dim,1:dim), &
                dim)
  end do

  do i=1,dim
    do j=1,dim
      write(6,'(a,i2,a,i2,a,2ES24.14)') 'SPIN-PHASE:  Zout(',i,',',j,') = ',Zout(i,j)
    end do
  end do

  write(6,'(//)')
  do i=1,dim
    do j=1,dim
      if ((j == i-1) .or. (j == i+1)) then
        write(6,'(A,i2,A,i2,A, 3(2ES24.14,3x))') 'SPIN-PHASE: PHS(',i,',',j,') = (x,y,z) =',(phs(l,i,j),l=1,3)
      else
        cycle
      end if
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
