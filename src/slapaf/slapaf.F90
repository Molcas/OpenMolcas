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

subroutine SlapAf(ireturn)

use Slapaf_Info, only: CallLast, isFalcon
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "warnings.h"
integer(kind=iwp) :: Columbus, iStop
character(len=8) :: ELOOP

!                                                                      *
!***********************************************************************
!                                                                      *
! Let us go to work !

iStop = -1
call RlxCtl(iStop)
! in case it is dmrg calculation
!      call read_dmrg_parameter_for_mclr()
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue
!                                                                      *
!***********************************************************************
!                                                                      *
! Put out the proper return code.

if (iStop == 2) then
  !*********** columbus interface **************************************
  ! prevent slapaf to startup alaska
  ! The last MRCI or MRAQCC energy is written to the RUNFILE by the
  ! Columbus modules, so no need to recompute it internally

  call Get_iScalar('Columbus',Columbus)
  if (Columbus == 1) then
    iReturn = _RC_ALL_IS_WELL_
  else

    ! Compute the last energy (module Last_Energy)

    if (CallLast) then
      call Start_Last_Energy()

      iReturn = _RC_INVOKED_OTHER_MODULE_
    else
      iReturn = _RC_ALL_IS_WELL_
    end if
  end if

else if (iStop == 3) then

  ! Invoke Alaska

  call Start_Alaska()

  iReturn = _RC_INVOKED_OTHER_MODULE_

else if (iStop == 6) then

  ! Saddle TS optimization

  iReturn = _RC_INVOKED_OTHER_MODULE_

else if (iStop == 8) then

  ! Predicted energy change too large

  iReturn = _RC_NOT_CONVERGED_

else if (iStop == 16) then

  ! Optimization didn't converge!

  iReturn = _RC_NOT_CONVERGED_

else if (iStop == 1) then

  ! Continue looping!

  call GetEnvF('EMIL_InLoop',ELOOP)
  if (ELOOP == ' ') ELOOP = '0'
  if (ELOOP(1:1) /= '0') then
    iReturn = _RC_CONTINUE_LOOP_
  else
    iReturn = _RC_ALL_IS_WELL_
  end if
  if (isFalcon) iReturn = _RC_ALL_IS_WELL_

else if (iStop == 0) then

  ! ???

  iReturn = _RC_ALL_IS_WELL_

else

  ! Bad bad boy!

  iReturn = _RC_INTERNAL_ERROR_
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine SlapAf
