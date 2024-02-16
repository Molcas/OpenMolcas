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

use Definitions, only: u5, u6

implicit none
integer, intent(inout) :: nneq
integer, intent(inout) :: neq(nneq), nexch(nneq)
!local variables:
integer :: i, LineNr
character(len=72) :: LINE
logical :: ab_initio_all
logical :: DBG

DBG = .false.
if (DBG) write(u6,'(A)') 'Enter fetch_neq'
if (DBG) write(u6,'(A,i3)') 'fetch_neq:  nneq=',nneq

neq = 0
nexch = 0
!=========== End of default settings====================================
rewind(u5)
50 read(u5,'(A72)',end=998) LINE
if (DBG) write(u6,'(A)') LINE
call NORMAL(LINE)
if (LINE(1:5) /= '&POLY') Go To 50
LINENR = 0

100 call xFlush(u6)
read(u5,'(A72)',end=998) LINE
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') Go To 100
if (LINE == ' ') Go To 100
if (LINE(1:3) == 'END') Go To 210 !End

if (LINE(1:4) == 'NNEQ') then
  ! number of non-equivalent centers; type of all centers
  read(u5,*,err=997) NNEQ,ab_initio_all
  if (DBG) write(u6,'(A,i4,A,L2)') 'NNEQ=',NNEQ,' ab_initio_all=',ab_initio_all
  ! number of equivalent centers of type "i"
  read(u5,*,err=997) (NEQ(i),i=1,Nneq)
  if (DBG) write(u6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
  ! number of RASSI wf for exchange
  read(u5,*,err=997) (Nexch(i),i=1,Nneq)
  if (DBG) write(u6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)
  Go To 210
end if
Go To 100

210 continue

goto 999 ! end

!-----------------------------------------------------------------------
997 continue
write(u6,*) ' READIN: Error reading "poly_aniso.input" '
write(u6,*) ' near line nr.',LINENR+1
Go To 999
998 continue
write(u6,*) ' READIN: Unexpected End of input file.'
!-----------------------------------------------------------------------
999 continue
if (DBG) write(u6,'(A)') 'Exit fetch_neq'

return

end subroutine fetch_neq
