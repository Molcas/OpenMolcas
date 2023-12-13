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

module abdata

use Definitions, only: wp

implicit none
private

real(kind=wp), allocatable :: atab(:,:), btab(:,:), p0(:), tvalue(:)

public :: atab, btab, close_abdata, p0, read_abdata, tvalue

contains

subroutine read_abdata()

  use Definitions, only: iwp
  use stdalloc, only: mma_allocate

  integer(kind=iwp) :: i, k, lu_abdata, maxdeg, ntab1, ntab2
  character(len=8) :: key
  logical(kind=iwp) :: found_abdata
  character(len=*), parameter :: ABDATA_NAME = 'ABDATA'
  integer(kind=iwp), external :: isFreeUnit

# include "macros.fh"

  call f_Inquire(ABDATA_NAME,found_abdata)
  if (.not. found_abdata) then
    call warningmessage(2,' the abdata file does not exist.')
    call abend()
  end if
  lu_abdata = isFreeUnit(22)
  call molcas_open(lu_abdata,ABDATA_NAME)

  do
    read(lu_abdata,'(a8)') key
    if (key == 'NTAB1, N') exit
  end do
  read(lu_abdata,*) ntab1,ntab2,maxdeg
  call mma_allocate(atab,[0,maxdeg],[ntab1,ntab2],label='atab')
  call mma_allocate(btab,[0,maxdeg],[ntab1,ntab2],label='btab')
  call mma_allocate(p0,[ntab1,ntab2],label='p0')
  call mma_allocate(tvalue,[ntab1,ntab2],label='tvalue')
  do i=ntab1,ntab2
    do
      read(lu_abdata,'(a8)') key
      if (key == 'TAB POIN') exit
    end do
    read(lu_abdata,*) k,tvalue(i),p0(i)
    unused_var(k)
    read(lu_abdata,*)
    read(lu_abdata,*) atab(:,i)
    read(lu_abdata,*)
    read(lu_abdata,*) btab(:,i)
  end do

  close(lu_abdata)

  return

end subroutine read_abdata

subroutine close_abdata()

  use stdalloc, only: mma_deallocate

  if (allocated(atab)) call mma_deallocate(atab)
  if (allocated(btab)) call mma_deallocate(btab)
  if (allocated(p0)) call mma_deallocate(p0)
  if (allocated(tvalue)) call mma_deallocate(tvalue)

end subroutine close_abdata

end module abdata
