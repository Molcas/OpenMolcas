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

subroutine Cho_TranslateErrorCode(iErr,MolcasCode)

implicit none
integer iErr, MolcasCode
#include "warnings.h"

if (iErr == 3) then ! special code used in parallel
  MolcasCode = _RC_NOT_AVAILABLE_
else if (iErr == 100) then
  MolcasCode = _RC_CHO_DUM_
else if (iErr == 101) then
  MolcasCode = _RC_CHO_MEM_
else if (iErr == 102) then
  MolcasCode = _RC_CHO_INI_
else if (iErr == 103) then
  MolcasCode = _RC_CHO_LOG_
else if (iErr == 104) then
  MolcasCode = _RC_CHO_RUN_
else if (iErr == 105) then
  MolcasCode = _RC_CHO_INP_
else
  MolcasCode = _RC_GENERAL_ERROR_
end if

end subroutine Cho_TranslateErrorCode
