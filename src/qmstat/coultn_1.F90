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

! s-s interaction, normal case.
real*8 function CoulTN_1(RA,RB,C,dSepInv,ExpA,ExpB)

implicit real*8(a-h,o-z)

T1 = 0.25d0*(2.0d0+C)
T2 = 0.25d0*RA
TA = (1.0d0-C)**2*(T1+T2)*ExpA
T1 = 0.25d0*(2.0d0-C)
T2 = 0.25d0*RB
TB = (1.0d0+C)**2*(T1+T2)*ExpB
CoulTN_1 = dSepInv*(1.0d0-TA-TB)

return

end function CoulTN_1
