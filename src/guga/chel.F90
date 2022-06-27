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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine CHEL(IA,IB,IIM,IEL,ISTOP)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IA, IB, IIM, IEL
integer(kind=iwp), intent(out) :: ISTOP
integer(kind=iwp) :: IR, IRR

IR = IIM-1
! CHECK FOR A=0 , B=IEL
IRR = IR-IA
if (IRR >= 0) then
  if (IRR >= IB-IEL) then
    ISTOP = 0
    return
  end if
end if
if (IEL /= 1) then
  ! CHECK FOR A=1 , B=IEL-2
  IRR = IR-IA+1
  if (IRR >= 0) then
    if (IRR >= IB-IEL+2) then
      ISTOP = 0
      return
    end if
  end if
end if
ISTOP = 1

return

end subroutine CHEL
