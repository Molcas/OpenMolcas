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

subroutine pseudospin(M,d,z,iDir,iOpt,iprint)
! d - dimension (input)
! moment(l,d,d) (input)
! z - pseudospin eigenfunctions (output)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: d, iDir, iOpt, iprint
complex(kind=wp), intent(in) :: M(3,d,d)
complex(kind=wp), intent(out) :: z(d,d)
integer(kind=iwp) :: i, info
real(kind=wp), allocatable :: w(:)
complex(kind=wp), allocatable :: M_tmp(:,:), z1(:,:)
real(kind=wp), external :: dznrm2_
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

call mma_allocate(W,d,'W')
call mma_allocate(Z1,d,d,'Z1')
call mma_allocate(M_tmp,d,d,'M_tmp')
info = 0
M_tmp(:,:) = M(iDir,:,:)
call diag_c2(M_tmp,d,info,w,z1)
if (iprint >= 3) then
  do i=1,d
    write(u6,'(A,i3,A,F24.14)') 'i=',i,' eigenvalue=',w(i)
  end do
end if
if (info /= 0) then
  write(u6,'(5x,a)') 'PSEUDO::  diagonalization of the zeeman hamiltonian failed.'
else
  if (iprint >= 3) then
    write(u6,*) 'PSEUDO:  norm of  M is:',dznrm2_(3*d*d,M,1)
    write(u6,*) 'PSEUDO:  norm of Z1 is:',dznrm2_(d*d,Z1,1)
  end if
  if (iDir == 3) then
    if (iOpt == 1) then
      call spin_phase(M,d,z1,z)
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
call mma_deallocate(M_tmp)

return

end subroutine pseudospin
