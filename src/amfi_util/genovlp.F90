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

subroutine genovlp(Lhigh,coulovlp,eval)
!bs generates overlap of normalized  primitives.

use AMFI_global, only: Lmax, MxprimL, normovlp, nprimit, OVLPinv, rootOVLP, rootOVLPinv
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lhigh
real(kind=wp), intent(in) :: coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)
real(kind=wp), intent(out) :: eval(MxprimL)
integer(kind=iwp) :: ipnt, Irun, Jrun, L, n
real(kind=wp) :: fact
real(kind=wp), allocatable :: evecinv(:,:)
real(kind=wp), allocatable, target :: scratch(:)
real(kind=wp), pointer :: tmp(:,:)

call mma_allocate(evecinv,MxprimL,MxprimL,label='evecinv')
call mma_allocate(scratch,MxprimL**2,label='scratch')
tmp(1:MxprimL,1:MxprimL) => scratch

do L=0,Lhigh
  n = nprimit(L)
  normovlp(1:n,1:n,L) = coulovlp(1:n,1:n,0,0,L,L)
  !bs invert the matrix, not very elegant, but sufficient
  ipnt = 0
  do jrun=1,n
    do irun=1,jrun
      ipnt = ipnt+1
      scratch(ipnt) = normovlp(irun,jrun,L)
    end do
  end do
  evecinv(:,1:n) = Zero
  do Jrun=1,n
    evecinv(jrun,jrun) = One
  end do
  call Jacob(scratch,evecinv,n,MxprimL)
  do irun=1,n
    eval(irun) = sqrt(scratch((irun*(irun+1))/2))
  end do
  !bs ensure normalization of the vectors.
  do IRUN=1,n
    fact = Zero
    do JRUN=1,n
      fact = fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
    end do
    fact = One/sqrt(fact)
    evecinv(:,IRUN) = fact*evecinv(:,IRUN)
  end do
  !bs now generate rootOVLP
  do irun=1,n
    tmp(1:n,irun) = eval(irun)*evecinv(1:n,irun)
  end do
  call dgemm_('N','T',n,n,n,One,evecinv,MxprimL,tmp,MxprimL,Zero,rootOVLP(:,:,L),MxprimL)
  !bs now generate rootOVLPinv
  do irun=1,n
    tmp(1:n,irun) = evecinv(1:n,irun)/eval(irun)
  end do
  call dgemm_('N','T',n,n,n,One,evecinv,MxprimL,tmp,MxprimL,Zero,rootOVLPinv(:,:,L),MxprimL)
  !bs now generate OVLPinv
  do irun=1,n
    tmp(1:n,irun) = evecinv(1:n,irun)/eval(irun)**2
  end do
  call dgemm_('N','T',n,n,n,One,evecinv,MxprimL,tmp,MxprimL,Zero,OVLPinv(:,:,L),MxprimL)
end do

nullify(tmp)
call mma_deallocate(evecinv)
call mma_deallocate(scratch)

return

end subroutine genovlp
