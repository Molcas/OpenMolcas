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

subroutine ChkInp_ccsort()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Check input for consistency                                      *
!                                                                      *
!***********************************************************************

use ccsort_global, only: IPT2, NBAS, nBasX, NCONF, NSYM, nSymX
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iErr, iSym

! Just print warning...
if (IPT2 == 0) then
  write(u6,*)
  write(u6,*) '       !!!!!WARNING!!!!!'
  write(u6,*)
  write(u6,*) '      *** input error ***'
  write(u6,*) '  The JOBIPH file does not include canonical orbitals'
  write(u6,*)
  write(u6,*) '       !!!!!WARNING!!!!!'
  write(u6,*)
  !call Quit_OnUserError()
end if

if (NCONF /= 1) then
  write(u6,*)
  write(u6,*) '  *** input error ***'
  write(u6,*) '  The JOBIPH file does not include a RHF or ROHF wave function'
  write(u6,*)
  call Quit_OnUserError()
end if

iErr = 0
if (nSym /= nSymX) iErr = 1
do iSym=1,nSym
  if (nBas(iSym) /= nBasX(iSym)) then
    iErr = 1
    exit
  end if
end do
if (iErr /= 0) then
  write(u6,*)
  write(u6,*) '  *** input error ***'
  write(u6,*) '  The JOBIPH and the TRAONE files are inconsistent'
  write(u6,*)
  call Quit_OnUserError()
end if

return

end subroutine ChkInp_ccsort
