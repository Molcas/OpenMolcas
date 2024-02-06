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

subroutine read_aniso_old_exch(input_file_name,nss,eso,MM,MS,ML)

use Constants, only: Zero, cZero, gElectron
use Definitions, only: wp

implicit none
#include "stdalloc.fh"
integer, intent(in) :: nss
real(kind=8), intent(out) :: eso(nss)
complex(kind=8), intent(out) :: MM(3,nss,nss)
complex(kind=8), intent(out) :: MS(3,nss,nss)
complex(kind=8), intent(out) :: ML(3,nss,nss)
character(len=180) :: input_file_name
! local variables:
integer :: nss_local, nstate_local
integer :: l, j, j1, j2, LuAniso, IsFreeUnit
real(kind=8), parameter :: g_e = 2.00231930437180_wp !IFG -gElectron
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:), tmp(:)
external :: IsFreeUnit
! in this subroutine nss is input data

! set to zero all arrays:
eso(:) = Zero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
ML(:,:,:) = cZero
nss_local = 0
nstate_local = 0

! read the file "aniso.input":
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,trim(input_file_name))
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate_local,nss_local
!-----------------------------------------------------------------------
call mma_allocate(tmp,nss_local,'tmp')
call dcopy_(nss_local,[Zero],0,tmp,1)
! local spin-orbit energy
read(LuAniso,*) (tmp(j),j=1,nss_local)
! copy the lowest nss states to eso:
do j=1,nss
  eso(j) = tmp(j)
end do
call mma_deallocate(tmp)
!-----------------------------------------------------------------------
read(LuAniso,*) (l,j=1,nstate_local)
!-----------------------------------------------------------------------

call mma_allocate(tmpR,nss_local,nss_local,'tmpR')
call mma_allocate(tmpI,nss_local,nss_local,'tmpI')
! magnetic moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  do j1=1,nss_local
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss_local)
  end do
  do j1=1,nss
    do j2=1,nss
      MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do

! spin moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  do j1=1,nss_local
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss_local)
  end do
  do j1=1,nss
    do j2=1,nss
      MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

! compute the orbital moment
do l=1,3
  do j1=1,nss
    do j2=1,nss
      ML(l,j1,j2) = -MM(l,j1,j2)-MS(l,j1,j2)*g_e
    end do
  end do
end do

close(LuAniso)

return

end subroutine read_aniso_old_exch
