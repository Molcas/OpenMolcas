!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Ignacio Fdez. Galvan                             *
!***********************************************************************

module text_file

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
private

public :: extend_line, next_non_comment

contains

!=======================================================================
! Read the next non-comment line, no matter its length.
! Return .false. if the file ends.
function next_non_comment(lu,line)

logical(kind=iwp) :: next_non_comment
integer(kind=iwp), intent(in) :: lu
character(len=:), allocatable, intent(inout) :: line
integer(kind=iwp) :: stat

next_non_comment = .false.
do
  call next_line(lu,line,stat)
  if (is_iostat_end(stat)) then
    ! if file ended, return false
    exit
  else if (stat == 0) then
    ! if line is not empty or comment, return true
    ! if it is empty or comment, read another line
    line(:) = adjustl(line)
    if ((line(1:1) /= '*') .and. (line /= '')) then
      next_non_comment = .true.
      exit
    end if
  else
    ! any other error, abort
    call abend()
  end if
end do

return

end function next_non_comment

!=======================================================================
! Read the next line, no matter its length.
subroutine next_line(lu,line,stat)

integer(kind=iwp), intent(in) :: lu
character(len=:), allocatable, intent(inout) :: line
integer(kind=iwp), intent(out) :: stat
integer(kind=iwp) :: readl
character(len=128) :: buf

if (allocated(line)) call mma_deallocate(line)
do
  read(lu,'(A)',iostat=stat,advance='no',size=readl) buf
  if (is_iostat_eor(stat)) then
    ! if line ended, add the chunk and return success
    call extend_line(line,buf(:readl),chop=.true.)
    stat = 0
    exit
  else if (stat == 0) then
    ! if line didn't end, add the chunk and read next
    call extend_line(line,buf,chop=.true.)
  else
    ! any other error, return
    exit
  end if
end do

end subroutine next_line

!=======================================================================
! Concatenate two strings into a longer one
! (if chop, the last character of the first string is removed)
subroutine extend_line(dynline,line,chop)

character(len=:), allocatable, intent(inout) :: dynline
character(len=*), intent(in) :: line
logical(kind=iwp), optional, intent(in) :: chop
logical(kind=iwp) :: chop_
character(len=:), allocatable :: aux

chop_ = .false.
if (present(chop)) chop_ = chop

! adding always a space to avoid 0-length strings
if (allocated(dynline)) then
  if (chop_) then
    call mma_allocate(aux,len(dynline)+len(line),label='AuxLine')
    aux(:) = dynline(:len(dynline)-1)//line//' '
  else
    call mma_allocate(aux,len(dynline)+len(line)+1,label='AuxLine')
    aux(:) = dynline//line//' '
  end if
  call mma_deallocate(dynline)
  ! move_alloc does not work properly in all compilers
# ifdef _USE_MOVE_ALLOC_
  call move_alloc(aux,dynline)
# else
  call mma_allocate(dynline,len(aux),label='ExtLine')
  dynline(:) = aux
  call mma_deallocate(aux)
# endif
else
  call mma_allocate(dynline,len(line)+1,label='Ext2Line')
  dynline(:) = line//' '
end if

end subroutine extend_line

end module text_file
