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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ThisIsRestrictedCode(developer_name,code_name,Abort)
! Thomas Bondo Pedersen, December 2012.
!
! An extension of OnlyIMayUseIt, this routine allows you to
! identify the code portion which is not for others to use.
! Based on (and using) OnlyIMayUseIt by V. Veryazov.
!
! If Abort: stop the execution

implicit none
character*(*) developer_name
character*(*) code_name
logical Abort
#include "warnings.fh"

character*256 val

val = ' '
call GetEnvf('MOLCAS_ISDEV',val)
if (val == 'PRODUCTION') return

#if defined (_DEBUGPRINT_)
call OnlyIMayUseIt(developer_name)
write(6,'(A,A)') '>>>>> Restricted code: ',code_name
write(6,'(A,A,//)') '>>>>> Contact ',developer_name
if (Abort) then
  if (val == ' ' .or. val /= developer_name) then
    call xQuit(_RC_GENERAL_ERROR_)
  end if
end if
call xFlush(6)
#else
if (val == ' ' .or. val /= developer_name) then
  call OnlyIMayUseIt(developer_name)
  write(6,'(A,A,//)') '>>>>> Restricted code: ',code_name
  if (Abort) call xQuit(_RC_GENERAL_ERROR_)
  call xFlush(6)
end if
#endif

end subroutine ThisIsRestrictedCode
