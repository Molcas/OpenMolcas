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

subroutine read_formatted_aniso_poly_NEW(input_file_name,nss,nstate,eso,MM,MS,iReturn)

use Constants, only: Zero, cZero, auTocm
use Definitions, only: wp, u6

implicit none
#include "stdalloc.fh"
integer, intent(inout) :: nss, nstate, iReturn
! nLoc is the maximal value of the array nss(1:nneq)
real(kind=8), intent(out) :: eso(nss)
complex(kind=8), intent(out) :: MM(3,nss,nss)
complex(kind=8), intent(out) :: MS(3,nss,nss)
character(len=180) :: input_file_name
! local variables:
real(kind=wp), allocatable :: eso_au(:)
integer :: i, j, l, LuAniso, IsFreeUnit
external :: IsFreeUnit
logical :: dbg

dbg = .false.
iReturn = 0
eso(:) = Zero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
call mma_allocate(eso_au,nss,'eso_au')
eso_au(:) = Zero

LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,input_file_name)

call read_nss(LuAniso,nss,dbg)
if (dbg) write(u6,*) 'read_formatted_aniso_poly_NEW: nss=',nss
call read_nstate(LuAniso,nstate,dbg)
if (dbg) write(u6,*) 'read_formatted_aniso_poly_NEW: nstate=',nstate
call read_eso(LuAniso,nss,eso_au,dbg)
if (dbg) write(u6,*) 'read_formatted_aniso_poly_NEW: eso_au=',(eso_au(i),i=1,nss)
call read_magnetic_moment(LuAniso,nss,MM(1:3,1:nss,1:nss),dbg)
if (dbg) write(u6,*) 'Call read_spin_moment'
flush(u6)
call read_spin_moment(LuAniso,nss,MS,dbg)

! compute the relative spin-orbit energies in cm-1
do i=1,nss
  eso(i) = (eso_au(i)-eso_au(1))*auTocm
end do
call mma_deallocate(eso_au)
close(LuAniso)

if (dbg) then
  write(u6,*) 'read_formatted_aniso_poly_NEW:  nss: ',nss
  write(u6,*) 'read_formatted_aniso_poly_NEW:   MM: '
  do l=1,3
    write(u6,'(A,I0)') 'projection: L=',l
    do i=1,nss
      write(u6,'(10(2F8.4,2x))') (MM(l,i,j),j=1,nss)
    end do
  end do

  write(u6,*) 'read_formatted_aniso_poly_NEW:   MS'
  do l=1,3
    write(u6,'(A,I0)') 'projection: L=',l
    do i=1,nss
      write(u6,'(10(2F8.4,2x))') (MS(l,i,j),j=1,nss)
    end do
  end do
end if

return

end subroutine read_formatted_aniso_poly_NEW
