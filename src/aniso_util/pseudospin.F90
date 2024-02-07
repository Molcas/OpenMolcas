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

subroutine pseudospin(M,dim,z,iDir,iOpt,iprint)
! dim - dimension (input)
! moment(l,dim,dim) (input)
! z - pseuDospin eigenfunctions (output)

use Constants, only: Zero, cZero
use Definitions, only: u6

implicit none
#include "stdalloc.fh"
integer, intent(in) :: dim, iprint
integer, intent(in) :: iDir, iOpt
complex(kind=8), intent(in) :: M(3,dim,dim)
complex(kind=8), intent(out) :: z(dim,dim)
! local variables:
integer :: info, i
real(kind=8), allocatable :: w(:)
complex(kind=8), allocatable :: z1(:,:)
real(kind=8) :: dznrm2_
external :: dznrm2_
logical :: dbg
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

call mma_allocate(W,dim,'W')
call mma_allocate(Z1,dim,dim,'Z1')
dbg = iprint >= 3
W(:) = Zero
Z(:,:) = cZero
Z1(:,:) = cZero
info = 0
call diag_c2(M(iDir,:,:),dim,info,w,z1)
if (dbg) then
  do i=1,dim
    write(u6,'(A,i3,A,F24.14)') 'i=',i,' eigenvalue=',w(i)
  end do
end if
if (info /= 0) then
  write(u6,'(5x,a)') 'PSEUDO::  diagonalization of the zeeman hamiltonian failed.'
else
  if (dbg) then
    write(u6,*) 'PSEUDO:  norm of  M is:',dznrm2_(3*dim*dim,M,1)
    write(u6,*) 'PSEUDO:  norm of Z1 is:',dznrm2_(dim*dim,Z1,1)
  end if
  if (iDir == 3) then
    if (iOpt == 1) then
      call spin_phase(M,dim,z1,z)
    else
      z(:,:) = z1(:,:)
      write(u6,*) 'PSEUDOSPIN:  iOpt = ',iOpt
      call WarningMessage(2,'PSEUDOSPIN: iOpt is not understood.')
    end if
  else
    z(:,:) = z1(:,:)
  end if
end if

call mma_deallocate(W)
call mma_deallocate(Z1)

return

end subroutine pseudospin
