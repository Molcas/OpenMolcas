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

subroutine fetch_data_RunFile_init(nss,nstate)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nss, nstate
integer(kind=iwp) :: mxjob, ndata, njob
logical(kind=iwp) :: FOUND

! check the presence of the RUNFILE
call f_inquire('RUNFILE',FOUND)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The RUNFILE was not found in the $WorkDir'
  write(u6,'(5X,A)') 'Are you running this calculation in the'
  write(u6,'(5X,A)') 'same $WorkDir where RASSI calculation'
  write(u6,'(5X,A)') 'was executed?'
  write(u6,'(5X,A)') 'Check your calculation again, and If'
  write(u6,'(5X,A)') 'necessary, submit a BUG report.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

! check the presence of the necessary information on RUNFILE:
call qpg_iscalar('NSS_SINGLE',FOUND)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The NSS Value was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_iscalar('NJOB_SINGLE',FOUND)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The NJOB Value was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_iscalar('MXJOB_SINGLE',FOUND)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The MXJOB Value was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_iscalar('NSTATE_SINGLE',FOUND)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The NSTATE Value was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

! fetch scalar data from RUNFILE
call get_iscalar('NSS_SINGLE',NSS)
call get_iscalar('NJOB_SINGLE',NJOB)
call get_iscalar('MXJOB_SINGLE',MXJOB)
call get_iscalar('NSTATE_SINGLE',NSTATE)

! check the presence of saved arrays on RUNFILE:

call qpg_iArray('MLTP_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The MLTP array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_iArray('JBNUM_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The JBNUM array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_iArray('LROOT_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The LROOT array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_dArray('ESO_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The ESO array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_dArray('UMATR_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The UMATR array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_dArray('UMATI_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The UMATI array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Please report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_dArray('ANGM_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The ANGMOM array was not found on RUNFILE'
  write(u6,'(5X,A)') '1. Check If ANGM keyword was used for SEWARD.'
  write(u6,'(5X,A)') '2. Check If MEES keyword was used for RASSI.'
  write(u6,'(5X,A)') '3. Check If PROP keyword was used for RASSI:'
  write(u6,'(9X,A)') 'PROP'
  write(u6,'(9X,A)') '3'
  write(u6,'(9X,A)') ' ''ANGMOM'' 1'
  write(u6,'(9X,A)') ' ''ANGMOM'' 2'
  write(u6,'(9X,A)') ' ''ANGMOM'' 3'
  write(u6,'(5X,A)') 'If MEES, ANGMOM and PROP keywords  were used and you still see this problem,'
  write(u6,'(5X,A)') 'please, report a BUG.'
  call xFlush(u6)
  call Quit_OnUserError()
end if

call qpg_dArray('DIP1_SINGLE',FOUND,NDATA)
if (.not. FOUND) then
  write(u6,'(5X,A)') 'The DIPMOM array was not found on RUNFILE'
  write(u6,'(5X,A)') 'Absorption intensities will not be computed'
end if

return

end subroutine fetch_data_RunFile_init
