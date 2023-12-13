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

subroutine chpcmp2_cvb(itst,iret)

use casvb_global, only: iprm, lstprm, mxprm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: itst
integer(kind=iwp), intent(out) :: iret

iprm = iprm+1
if (iprm > mxprm) then
  write(u6,*) ' Dimensioning error in CHPCMP2!',iprm,mxprm
  call abend_cvb()
end if
iret = lstprm(iprm)
lstprm(iprm) = itst

return

end subroutine chpcmp2_cvb
