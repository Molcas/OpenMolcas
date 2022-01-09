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


subroutine LDF_SetOptionFlag(Option,value)
implicit none
character*4 Option
logical value
#include "localdf.fh"

if (Option == 'LDF2') then
  LDF2 = value
else if (Option == 'CHEC') then
  CheckPairIntegrals = value
else if (Option == 'VERI') then
  VerifyFit = value
else if (Option == 'OVER') then
  CheckOverlapIntegrals = value
else if (Option == 'WRUC') then
  WriteUnconstrainedC = value
else if (Option == 'UNIQ') then
  UseUniqueAtomPairs = value
else
  call WarningMessage(2,'LDF_SetOptionFlag: unknown Option')
  write(6,'(A,A)') 'Option=',Option
  write(6,'(A,L1)') 'Value=',value
  call LDF_Quit(1)
end if

return

end subroutine LDF_SetOptionFlag
