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

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(inout) :: nss, nstate
integer, intent(out) :: multiplicity(nstate)
real(kind=8), intent(out) :: eso(nss)
complex(kind=8), intent(out) :: MM(3,nss,nss)
complex(kind=8), intent(out) :: MS(3,nss,nss)
complex(kind=8), intent(out) :: ML(3,nss,nss)
character(Len=180) :: input_file_name
! local variables:
integer :: l, j, j1, j2, LuAniso, IsFreeUnit
real(kind=8) :: g_e
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
external :: IsFreeUnit

g_e = 2.00231930437180_wp
!   set to zero all arrays:
multiplicity = 0
eso = 0.0_wp
MM = (0.0_wp,0.0_wp)
MS = (0.0_wp,0.0_wp)
ML = (0.0_wp,0.0_wp)
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
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do

! spin moment
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
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

!if (iprint > 4) then
!  write(6,'(10a12)') (('------------'),j=1,10)
!  write(6,'(15x,a,i2,a)') 'aniso_.input'
!  write(6,'(10a12)') (('------------'),j=1,10)
!  write(6,'(5x,a,i6)') 'nstate = ',nstate
!  write(6,'(5x,a,i6)') '   nss = ',nss
!  write(6,'(5x,a)') ' eso(j): '
!  write(6,'(10(f12.5,1x))') (eso(j),j=1,nss)
!  write(6,'(5x,a,i2,a)') 'multiplicity(j):'
!  write(6,'(40i3)') (multiplicity(j),j=1,nstate)
!  write(6,'(5x,a,i2,a)') 'dipso(l,j1,j2):'
!  do l=1,3
!    write(6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(6,'(20(2f20.14,1x))') (MM(l,j1,j2),j2=1,nss)
!    end do
!  end do
!  write(6,*)
!  write(6,'(5x,a,i2,a)') 's_so(l,j1,j2):'
!  do l=1,3
!    write(6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(6,'(20(2f20.14,1x))') (MS(l,j1,j2),j2=1,nss)
!    end do
!  end do
!end if
close(LuAniso)

return

end subroutine read_formatted_aniso_old
