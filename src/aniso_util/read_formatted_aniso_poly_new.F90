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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, auTocm
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(inout) :: nss, nstate
real(kind=wp), intent(out) :: eso(nss)
complex(kind=wp), intent(out) :: MM(3,nss,nss), MS(3,nss,nss)
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: LuAniso
real(kind=wp), allocatable :: eso_au(:)
complex(kind=wp), allocatable :: MTMP(:,:,:)
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
integer(kind=iwp) :: i, j, l
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_
integer(kind=iwp), external :: IsFreeUnit

iReturn = 0
MM(:,:,:) = cZero
MS(:,:,:) = cZero
call mma_allocate(eso_au,nss,'eso_au')
eso_au(:) = Zero

LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,input_file_name)

call read_nss(LuAniso,nss,dbg)
#ifdef _DEBUGPRINT_
write(u6,*) 'read_formatted_aniso_poly_NEW: nss=',nss
#endif
call read_nstate(LuAniso,nstate,dbg)
#ifdef _DEBUGPRINT_
write(u6,*) 'read_formatted_aniso_poly_NEW: nstate=',nstate
#endif
call read_eso(LuAniso,nss,eso_au,dbg)
#ifdef _DEBUGPRINT_
write(u6,*) 'read_formatted_aniso_poly_NEW: eso_au=',(eso_au(i),i=1,nss)
#endif
call mma_allocate(MTMP,3,nss,nss)
call read_magnetic_moment(LuAniso,nss,MTMP,dbg)
MM(:,1:nss,1:nss) = MTMP(:,:,:)
#ifdef _DEBUGPRINT_
write(u6,*) 'Call read_spin_moment'
call xFlush(u6)
#endif
call read_spin_moment(LuAniso,nss,MTMP,dbg)
MS(:,1:nss,1:nss) = MTMP(:,:,:)
call mma_deallocate(MTMP)

! compute the relative spin-orbit energies in cm-1
eso(:) = (eso_au(:)-eso_au(1))*auTocm
call mma_deallocate(eso_au)
close(LuAniso)

#ifdef _DEBUGPRINT_
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
#endif

return

end subroutine read_formatted_aniso_poly_NEW
