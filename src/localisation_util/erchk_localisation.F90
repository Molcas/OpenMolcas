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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ERChk_Localisation(irc,lnBas,lnOcc,lnFro,lnSym)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lnSym, lnBas(lnSym), lnOcc(lnSym), lnFro(lnSym)
#include "cholesky.fh"
#include "choorb.fh"
integer(kind=iwp) :: iSym, nTst

irc = 0

if ((lnSym < 1) .or. (lnSym > 8)) then
  irc = 1
  return
end if

if (lnSym /= nSym) then
  irc = 2
  return
end if

do iSym=1,nSym
  if (lnBas(iSym) /= nBas(iSym)) then
    irc = 3
    return
  end if
  nTst = lnOcc(iSym)+lnFro(iSym)
  if (nTst > nBas(iSym)) then
    irc = 4
    return
  end if
end do

end subroutine ERChk_Localisation
