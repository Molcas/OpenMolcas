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
! Copyright (C) 1996, Bernd Schimmelpfennig                            *
!***********************************************************************

function LMdepang(L,M,l1,l2,l3,l4,m1,m2,m3,m4,cheater)
!bs   l1-l4 and m1-m4 are already shifted !!
!bs   purpose: calculates the angular part of the
!bs   coulomb-type integrals. See documentation for details...
!bs   LMdepang= LM dependent angular factors
!bs   cheater included for a correcting signs, as there were some
!bs   signs (only signs!!!!) missing when compared to HERMIT
!bs                                        B.S.  08.10.96

use Constants, only: Zero, One, Four, Pi
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: LMdepang
integer(kind=iwp), intent(in) :: L, M, l1, l2, l3, l4, m1, m2, m3, m4
real(kind=wp), intent(in) :: cheater
integer(kind=iwp) :: isum
real(kind=wp) :: fact1, fact2, fact3, sgn
real(kind=wp), external :: couple3J

LMdepang = Zero
!bs some quick checks
if (L < abs(M)) return
if (l1 < abs(m1)) return
if (l2 < abs(m2)) return
if (l3 < abs(m3)) return
if (l4 < abs(m4)) return
!bs prefactor
fact1 = Four*Pi/real(L+L+1,wp)
!bs determining the sign
isum = -l3-l1-l4-l2+2*(M+m3+m4)   !???? I am not sure
if (mod(isum,4) == 0) then
  sgn = One
else if (abs(mod(isum,4)) == 2) then
  sgn = -One
else
  sgn = Zero
  write(u6,*) 'L,l1,l2,l3,l4,M,m1,m2,m3,m4'
  write(u6,'(10I3)') L,l1,l2,l3,l4,M,m1,m2,m3,m4
  write(u6,*) 'isum= ',isum,' mod = ',mod(isum,4)
  call SysHalt('lmdepang')
end if
fact2 = couple3J(L,l3,l1,-M,m3,-m1)
fact3 = couple3J(L,l4,l2,M,m4,-m2)
!write(6,*) 'fact2,fact3 ',fact2,fact3
LMdepang = cheater*sgn*fact1*fact2*fact3

return

end function LMdepang
