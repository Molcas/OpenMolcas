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
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Print_Isotopes()

use Basis_Info, only: dbsc, nCnttp
use Constants, only: UtoAU
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: i, iAtom, iPrint, iRout
real(kind=wp) :: act_Mass, def_Mass
logical(kind=iwp) :: Changed
real(kind=wp), external :: rMass

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Determine whether any mass is different from default

Changed = .false.
do i=1,nCnttp
  if (dbsc(i)%Aux .or. dbsc(i)%Frag) cycle
  iAtom = dbsc(i)%AtmNr
  if (dbsc(i)%CntMass /= rMass(iAtom)) then
    Changed = .true.
    exit
  end if
end do
if ((.not. Changed) .and. (iPrint <= 5)) return
!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,*)
call CollapseOutput(1,'   Isotope specification:')
write(u6,'(3X,A)') '   ----------------------'
write(u6,*)
if (Changed) then
  write(u6,10) 'Center                     [     Default     ]'
  write(u6,10) 'Type   Z    A    mass (Da) [   A    mass (Da)]'
  write(u6,10) '---------------------------------------------'
else
  write(u6,10) 'Center'
  write(u6,10) 'Type   Z    A    mass (Da)'
  write(u6,10) '--------------------------'
end if
do i=1,nCnttp
  if (dbsc(i)%Aux .or. dbsc(i)%Frag) cycle
  iAtom = dbsc(i)%AtmNr
  act_Mass = dbsc(i)%CntMass/UToAU
  def_Mass = rmass(iAtom)/UToAU
  if (act_Mass /= def_Mass) then
    write(u6,101) i,iAtom,nint(act_Mass),act_Mass,nint(def_Mass),def_Mass
  else
    write(u6,100) i,iAtom,nint(act_Mass),act_Mass
  end if
end do
call CollapseOutput(0,'   Isotope specification:')
write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *

return

10 format(1X,A)
100 format(I5,1X,I3,1X,I4,1X,F12.6)
101 format(I5,1X,I3,1X,I4,1X,F12.6,1X,'[',I4,1X,F12.6,']')

end subroutine Print_Isotopes
