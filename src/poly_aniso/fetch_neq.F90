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

subroutine fetch_neq(nneq,neq,nexch)

use Definitions, only: iwp, u5, u6

implicit none
integer(kind=iwp), intent(inout) :: nneq, neq(nneq), nexch(nneq)
integer(kind=iwp) :: i, istatus, LineNr
logical(kind=iwp) :: ab_initio_all
character(len=72) :: LINE

#include "macros.fh"

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'Enter fetch_neq'
write(u6,'(A,i3)') 'fetch_neq:  nneq=',nneq
#endif

neq = 0
nexch = 0
!=========== End of default settings====================================
rewind(u5)
do
  read(u5,'(A72)',iostat=istatus) LINE
  if (istatus < 0) then
    call Error(2)
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Exit fetch_neq'
#   endif
    return
  end if
# ifdef _DEBUGPRINT_
  write(u6,'(A)') LINE
# endif
  call NORMAL(LINE)
  if (LINE(1:5) == '&POLY') exit
end do
LINENR = 0

do
  call xFlush(u6)
  read(u5,'(A72)',iostat=istatus) LINE
  if (istatus < 0) then
    call Error(2)
    exit
  end if
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle
  if (LINE(1:3) == 'END') exit

  if (LINE(1:4) == 'NNEQ') then
    ! number of non-equivalent centers; type of all centers
    read(u5,*,iostat=istatus) NNEQ,ab_initio_all
    if (istatus /= 0) then
      call Error(1)
      exit
    end if
    unused_var(ab_initio_all)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,i4,A,L2)') 'NNEQ=',NNEQ,' ab_initio_all=',ab_initio_all
#   endif
    ! number of equivalent centers of type "i"
    read(u5,*,iostat=istatus) (NEQ(i),i=1,Nneq)
    if (istatus /= 0) then
      call Error(1)
      exit
    end if
#   ifdef _DEBUGPRINT_
    write(u6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
#   endif
    ! number of RASSI wf for exchange
    read(u5,*,iostat=istatus) (Nexch(i),i=1,Nneq)
    if (istatus /= 0) then
      call Error(1)
      exit
    end if
#   ifdef _DEBUGPRINT_
    write(u6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)
#   endif
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'Exit fetch_neq'
#endif

return

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      write(u6,*) ' FETCH_NEQ: Error reading "poly_aniso.input" '
      write(u6,*) ' near line nr.',LINENR+1
    case (2)
      write(u6,*) ' FETCH_NEQ: Unexpected End of input file.'
  end select

end subroutine Error

end subroutine fetch_neq
