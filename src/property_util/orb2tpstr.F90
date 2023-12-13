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

subroutine orb2tpstr(NSYM,NB,NF,NI,N1,N2,N3,NS,ND,TYPESTRING)

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSYM, NB(NSYM), NF(*), NI(*), N1(*), N2(*), N3(*), NS(*), ND(*)
character, intent(_OUT_) :: typestring(*)
integer(kind=iwp) :: iSym, iStart

iStart = 1
do isym=1,nsym
  call orb2tpstr_sym(NF(iSym),NI(iSym),N1(iSym),N2(iSym),N3(iSym),NS(iSym),ND(iSym),TYPESTRING(iStart))
  iStart = iStart+nB(iSym)
end do

end subroutine orb2tpstr
