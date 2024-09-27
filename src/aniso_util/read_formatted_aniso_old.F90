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

subroutine read_formatted_aniso_old(input_file_name,nss,nstate,multiplicity,eso,MM,MS,ML)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, gElectron
use Definitions, only: wp, iwp

implicit none
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(inout) :: nss, nstate
integer(kind=iwp), intent(out) :: multiplicity(nstate)
real(kind=wp), intent(out) :: eso(nss)
complex(kind=wp), intent(out) :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
integer(kind=iwp) :: j, j1, j2, l, LuAniso
real(kind=wp), allocatable :: tmpI(:,:), tmpR(:,:)
real(kind=wp), parameter :: g_e = -gElectron
integer(kind=iwp), external :: IsFreeUnit

! set to zero all arrays
multiplicity(:) = 0
eso(:) = Zero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
ML(:,:,:) = cZero
! read the file "aniso.input":
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,trim(input_file_name))
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate,nss
read(LuAniso,*) (eso(j),j=1,nss)
read(LuAniso,*) (multiplicity(j),j=1,nstate)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! magnetic moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! spin moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MS(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

! compute the orbital moment
ML(:,1:nss,1:nss) = -MM(:,1:nss,1:nss)-g_e*MS(:,1:nss,1:nss)

!if (iprint > 4) then
!  write(u6,'(10a12)') (('------------'),j=1,10)
!  write(u6,'(15x,a,i2,a)') 'aniso_.input'
!  write(u6,'(10a12)') (('------------'),j=1,10)
!  write(u6,'(5x,a,i6)') 'nstate = ',nstate
!  write(u6,'(5x,a,i6)') '   nss = ',nss
!  write(u6,'(5x,a)') ' eso(j): '
!  write(u6,'(10(f12.5,1x))') (eso(j),j=1,nss)
!  write(u6,'(5x,a,i2,a)') 'multiplicity(j):'
!  write(u6,'(40i3)') (multiplicity(j),j=1,nstate)
!  write(u6,'(5x,a,i2,a)') 'dipso(l,j1,j2):'
!  do l=1,3
!    write(u6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(u6,'(20(2f20.14,1x))') (MM(l,j1,j2),j2=1,nss)
!    end do
!  end do
!  write(u6,*)
!  write(u6,'(5x,a,i2,a)') 's_so(l,j1,j2):'
!  do l=1,3
!    write(u6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(u6,'(20(2f20.14,1x))') (MS(l,j1,j2),j2=1,nss)
!    end do
!  end do
!end if
close(LuAniso)

return

end subroutine read_formatted_aniso_old
