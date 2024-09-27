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
! in this subroutine nss is input data

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: gElectron
use Definitions, only: wp, iwp

implicit none
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(in) :: nss
real(kind=wp), intent(out) :: eso(nss)
complex(kind=wp), intent(out) :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
integer(kind=iwp) :: j, j1, j2, l, LuAniso, nss_local, nstate_local
real(kind=wp), allocatable :: tmp(:), tmpI(:,:), tmpR(:,:)
real(kind=wp), parameter :: g_e = -gElectron
integer(kind=iwp), external :: IsFreeUnit

! read the file "aniso.input":
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,trim(input_file_name))
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate_local,nss_local
!-----------------------------------------------------------------------
call mma_allocate(tmp,nss_local,'tmp')
! local spin-orbit energy
read(LuAniso,*) (tmp(j),j=1,nss_local)
! copy the lowest nss states to eso:
eso(:) = tmp(:)
call mma_deallocate(tmp)
!-----------------------------------------------------------------------
read(LuAniso,*) (l,j=1,nstate_local)
!-----------------------------------------------------------------------

call mma_allocate(tmpR,nss_local,nss_local,'tmpR')
call mma_allocate(tmpI,nss_local,nss_local,'tmpI')
! magnetic moment
do l=1,3
  do j1=1,nss_local
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss_local)
  end do
  MM(l,:,:) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! spin moment
do l=1,3
  do j1=1,nss_local
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss_local)
  end do
  MS(l,:,:) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

! compute the orbital moment
ML(:,:,:) = -MM(:,:,:)-g_e*MS(:,:,:)

close(LuAniso)

return

end subroutine read_aniso_old_exch
