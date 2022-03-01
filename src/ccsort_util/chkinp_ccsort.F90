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

implicit real*8(A-H,O-Z)
#include "ccsort.fh"
#include "motra.fh"

! Just print warning...
if (IPT2 == 0) then
  write(6,*)
  write(6,*) '       !!!!!WARNING!!!!!'
  write(6,*)
  write(6,*) '      *** input error ***'
  write(6,*) '  The JOBIPH file does not include canonical orbitals'
  write(6,*)
  write(6,*) '       !!!!!WARNING!!!!!'
  write(6,*)
  !call Quit_OnUserError()
end if

if (NCONF /= 1) then
  write(6,*)
  write(6,*) '  *** input error ***'
  write(6,*) '  The JOBIPH file does not include a RHF or ROHF wave function'
  write(6,*)
  call Quit_OnUserError()
end if

iErr = 0
if (nSym /= nSymX) iErr = 1
do iSym=1,nSym
  if (nBas(iSym) /= nBasX(iSym)) iErr = 1
end do
if (iErr /= 0) then
  write(6,*)
  write(6,*) '  *** input error ***'
  write(6,*) '  The JOBIPH and the TRAONE files are inconsistent'
  write(6,*)
  call Quit_OnUserError()
end if

return

end subroutine ChkInp_ccsort
