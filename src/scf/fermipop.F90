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
real*8 function FermiPop(e,o,n,T,nEle,UHF_occ)
!----------------------------------------------------------------------*
! Dummy arguments:                                                     *
! e(*) -- Orbital energies, input.                                     *
! o(*) -- Orbital occupations, output.                                 *
! n    -- Number of orbitals, input.                                   *
! T    -- Temperature, input.                                          *
! nEle -- Number of electrons.                                         *
!----------------------------------------------------------------------*

use Constants, only: Zero, Half, One, Three, Ten

implicit none
integer n, nEle
real*8 e(n), o(n), T, UHF_occ
!----------------------------------------------------------------------*
! Local variables:                                                     *
!----------------------------------------------------------------------*
real*8 ef, beta, f, f_old, Step, ff
real*8 x0, x1, x2, y0, y2, z
integer i, iter
#ifdef _DEBUGPRINT_
real*8 y1
#endif

!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
ef = Zero
if (T <= Zero) then
  beta = 1.0d99
else
  beta = One/T
end if
!----------------------------------------------------------------------*
! Scan for Fermi level                                                 *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(6,'(a)') 'Scan for Fermi energy level'
write(6,'(a)') '       ef             y       '
write(6,'(a)') ' -------------- --------------'
#endif

f = -nEle
f_old = f
do i=1,n
  !write(6,'(A,G20.10)') 'e(i)=',e(i)
  z = beta*(e(i)-ef)
  z = min(z,30.d0)
  f = f+UHF_occ/(One+exp(z))
end do
if (f > Zero) then
  Step = -One
else
  Step = One
end if
Iter = 0
100 continue
Iter = Iter+1
if (Iter > 100000) goto 101
f_old = f
ef = ef+Step
!f = -nEle
ff = Zero
!vv overoptimization with Intel compiler
i = 1
300 continue
!do i=1,n
z = beta*(e(i)-ef)
z = min(z,30.d0)
ff = ff+1/(One+exp(z))
i = i+1
if (i <= n) goto 300
!end do
f = -nEle+ff*UHF_occ
#ifdef _DEBUGPRINT_
write(6,'(2G20.10)') ef,f
#endif
if (f*f_old > Zero) goto 100
101 continue
!----------------------------------------------------------------------*
! Refine with interval halving.                                        *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(6,'(a)') 'Refine Fermi level with interval halving'
write(6,'(a)') '       y0            y2             y1       '
write(6,'(a)') ' -------------- -------------- --------------'
y1 = f
#endif
x0 = ef-Step
x1 = ef
y0 = f_old
x2 = half*(x0+x1)
Iter = 0
200 continue
Iter = iter+1
if (Iter > 1000) goto 201
ef = x2
f = -nEle
do i=1,n
  z = beta*(e(i)-ef)
  z = min(z,three*ten)
  f = f+UHF_occ/(One+exp(z))
end do
y2 = f
#ifdef _DEBUGPRINT_
write(6,'(3f15.8)') y0,y2,y1
#endif
if (abs(y2) < 1.0d-9) goto 201
if (y0*y2 <= Zero) then
  x1 = x2
# ifdef _DEBUGPRINT_
  y1 = y2
# endif
else
  x0 = x2
  y0 = y2
end if
x2 = half*(x0+x1)
goto 200
201 continue

!----------------------------------------------------------------------*
! Populate occupation number vector.                                   *
!----------------------------------------------------------------------*
!write(6,*)
f = Zero
do i=1,n
  !write(6,'(2G20.10)') e(i),ef
  z = beta*(e(i)-ef)
  !write(6,*) 'z,Beta=',z,Beta
  z = min(z,Three*Ten)
  !write(6,*) 'z=',z
  o(i) = UHF_occ/(One+exp(z))
  !write(6,'(1G20.10)') o(i)
  f = f+o(i)
end do
f = nEle/f
!write(6,*)
!write(6,'(1G20.10)') f
!write(6,*)
do i=1,n
  o(i) = f*o(i)
  !write(6,'(1G20.10)') o(i)
end do
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

FermiPop = ef

return

end function FermiPop
