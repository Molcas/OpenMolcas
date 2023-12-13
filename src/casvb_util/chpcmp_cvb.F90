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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

function chpcmp_cvb(itst)

use casvb_global, only: iprm, lstprm, mxprm
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: chpcmp_cvb
integer(kind=iwp), intent(in) :: itst

iprm = iprm+1
if (iprm > mxprm) then
  write(u6,*) ' Dimensioning error in CHPCMP!',iprm,mxprm
  call abend_cvb()
end if
chpcmp_cvb = (lstprm(iprm) /= itst)
lstprm(iprm) = itst

return

end function chpcmp_cvb
