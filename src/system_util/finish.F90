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
! Copyright (C) 2001-2016, Valera Veryazov                             *
!***********************************************************************

subroutine finish(rc)
! Gracefully shuts down a program module.
! After everything is closed properly, xquit is
! called to do the actual termination.

use warnings, only: MaxWarnMess
use Symmetry_Info, only: Symmetry_Info_Free
use Isotopes, only: Free_Isotopes
#ifndef _HAVE_EXTRA_
use Prgm, only: prgmfree
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: rc
integer(kind=iwp) :: idum = 0
#include "WrkSpc.fh"

call Symmetry_Info_Free()
call Free_Isotopes()

call fin_run_use()
#ifndef _HAVE_EXTRA_
call prgmfree()
#endif

call GetMem('ip_iDum','Free','Inte',ip_iDummy,1)
call GetMem('ip_Dum','Free','Real',ip_Dummy,1)
call GetMem('Finish','List','Real',iDum,iDum)
call GetMem('Finish','Term','Real',iDum,iDum)

call StatusLine('Happy landing',' ')
if (MaxWarnMess > 1) then
  call WarningMessage(1,'There were warnings during the execution;Please, check the output with care!')
end if

#ifdef _HAVE_EXTRA_
call prgmfree()
#endif
call AixCheck()
call xml_close('module')

#ifdef _DELAYED_
call close_BLAS()
#endif

call xquit(rc)

end subroutine finish
