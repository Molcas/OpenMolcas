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

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iErr
integer(kind=iwp), intent(out) :: MolcasCode
#include "warnings.h"

select case (iErr)
  case (3) ! special code used in parallel
    MolcasCode = _RC_NOT_AVAILABLE_
  case (100)
    MolcasCode = _RC_CHO_DUM_
  case (101)
    MolcasCode = _RC_CHO_MEM_
  case (102)
    MolcasCode = _RC_CHO_INI_
  case (103)
    MolcasCode = _RC_CHO_LOG_
  case (104)
    MolcasCode = _RC_CHO_RUN_
  case (105)
    MolcasCode = _RC_CHO_INP_
  case default
    MolcasCode = _RC_GENERAL_ERROR_
end select

end subroutine Cho_TranslateErrorCode
