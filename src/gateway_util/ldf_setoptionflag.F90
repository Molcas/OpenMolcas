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

subroutine LDF_SetOptionFlag(Option,Val)

use Definitions, only: iwp, u6

implicit none
character(len=4), intent(in) :: Option
logical(kind=iwp), intent(in) :: Val
#include "localdf.fh"

if (Option == 'LDF2') then
  LDF2 = Val
else if (Option == 'CHEC') then
  CheckPairIntegrals = Val
else if (Option == 'VERI') then
  VerifyFit = Val
else if (Option == 'OVER') then
  CheckOverlapIntegrals = Val
else if (Option == 'WRUC') then
  WriteUnconstrainedC = Val
else if (Option == 'UNIQ') then
  UseUniqueAtomPairs = Val
else
  call WarningMessage(2,'LDF_SetOptionFlag: unknown Option')
  write(u6,'(A,A)') 'Option=',Option
  write(u6,'(A,L1)') 'Val=',Val
  call LDF_Quit(1)
end if

return

end subroutine LDF_SetOptionFlag
