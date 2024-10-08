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
! Copyright (C) 1995, Martin Schuetz                                   *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This function computes the Fermi energy level for a number of        *
! energy levels and populates them. Each level is populated with up    *
! to 2 electrons.                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: June 1999                                                   *
! History:                                                             *
!                                                                      *
!***********************************************************************

!#define _DEBUGPRINT_
function FermiPop(e,o,n,T,nEle,UHF_occ)
!----------------------------------------------------------------------*
! Dummy arguments:                                                     *
! e(*) -- Orbital energies, input.                                     *
! o(*) -- Orbital occupations, output.                                 *
! n    -- Number of orbitals, input.                                   *
! T    -- Temperature, input.                                          *
! nEle -- Number of electrons.                                         *
!----------------------------------------------------------------------*

use Constants, only: Zero, One, Three, Ten, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp) :: FermiPop
integer(kind=iwp), intent(in) :: n, nEle
real(kind=wp), intent(in) :: e(n), T, UHF_occ
real(kind=wp), intent(out) :: o(n)
integer(kind=iwp) :: i, iter
real(kind=wp) :: beta, ef, f, f_old, ff, Step, x0, x1, x2, y0, y2, z
#ifdef _DEBUGPRINT_
real(kind=wp) :: y1
#endif

!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
ef = Zero
if (T <= Zero) then
  beta = 1.0e99_wp
else
  beta = One/T
end if
!----------------------------------------------------------------------*
! Scan for Fermi level                                                 *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(u6,'(a)') 'Scan for Fermi energy level'
write(u6,'(a)') '       ef             y       '
write(u6,'(a)') ' -------------- --------------'
#endif

f = -nEle
f_old = f
do i=1,n
  !write(u6,'(A,G20.10)') 'e(i)=',e(i)
  z = min(beta*(e(i)-ef),30.0_wp)
  f = f+UHF_occ/(One+exp(z))
end do
if (f > Zero) then
  Step = -One
else
  Step = One
end if
Iter = 0
do
  Iter = Iter+1
  if (Iter > 100000) exit
  f_old = f
  ef = ef+Step
  !f = -nEle
  ff = Zero
  do i=1,n
    z = min(beta*(e(i)-ef),30.0_wp)
    ff = ff+1/(One+exp(z))
  end do
  f = -nEle+ff*UHF_occ
# ifdef _DEBUGPRINT_
  write(u6,'(2G20.10)') ef,f
# endif
  if (f*f_old <= Zero) exit
end do
!----------------------------------------------------------------------*
! Refine with interval halving.                                        *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(u6,'(a)') 'Refine Fermi level with interval halving'
write(u6,'(a)') '       y0            y2             y1       '
write(u6,'(a)') ' -------------- -------------- --------------'
y1 = f
#endif
x0 = ef-Step
x1 = ef
y0 = f_old
x2 = half*(x0+x1)
Iter = 0
do
  Iter = iter+1
  if (Iter > 1000) exit
  ef = x2
  f = -nEle
  do i=1,n
    z = min(beta*(e(i)-ef),three*ten)
    f = f+UHF_occ/(One+exp(z))
  end do
  y2 = f
# ifdef _DEBUGPRINT_
  write(u6,'(3f15.8)') y0,y2,y1
# endif
  if (abs(y2) < 1.0e-9_wp) exit
  if (y0*y2 <= Zero) then
    x1 = x2
#   ifdef _DEBUGPRINT_
    y1 = y2
#   endif
  else
    x0 = x2
    y0 = y2
  end if
  x2 = half*(x0+x1)
end do

!----------------------------------------------------------------------*
! Populate occupation number vector.                                   *
!----------------------------------------------------------------------*
!write(u6,*)
f = Zero
do i=1,n
  !write(u6,'(2G20.10)') e(i),ef
  z = beta*(e(i)-ef)
  !write(u6,*) 'z,Beta=',z,Beta
  z = min(z,Three*Ten)
  !write(u6,*) 'z=',z
  o(i) = UHF_occ/(One+exp(z))
  !write(u6,'(1G20.10)') o(i)
  f = f+o(i)
end do
f = nEle/f
!write(u6,*)
!write(u6,'(1G20.10)') f
!write(u6,*)
o(:) = f*o(:)
!do i=1,n
!  write(u6,'(1G20.10)') o(i)
!end do
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

FermiPop = ef

return

end function FermiPop
