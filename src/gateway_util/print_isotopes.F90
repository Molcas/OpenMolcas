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

use Period
use Basis_Info, only: nCnttp, dbsc

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "constants2.fh"
logical Changed

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
LuWr = 6
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
write(LuWr,*)
call CollapseOutput(1,'   Isotope specification:')
write(LuWr,'(3X,A)') '   ----------------------'
write(LuWr,*)
if (Changed) then
  write(LuWr,10) 'Center                     [     Default     ]'
  write(LuWr,10) 'Type   Z    A    mass (Da) [   A    mass (Da)]'
  write(LuWr,10) '---------------------------------------------'
else
  write(LuWr,10) 'Center'
  write(LuWr,10) 'Type   Z    A    mass (Da)'
  write(LuWr,10) '--------------------------'
end if
do i=1,nCnttp
  if (dbsc(i)%Aux .or. dbsc(i)%Frag) cycle
  iAtom = dbsc(i)%AtmNr
  act_Mass = dbsc(i)%CntMass/UToAU
  def_Mass = rmass(iAtom)/UToAU
  if (act_Mass /= def_Mass) then
    write(LuWr,101) i,iAtom,nint(act_Mass),act_Mass,nint(def_Mass),def_Mass
  else
    write(LuWr,100) i,iAtom,nint(act_Mass),act_Mass
  end if
end do
call CollapseOutput(0,'   Isotope specification:')
write(LuWr,*)
!                                                                      *
!***********************************************************************
!                                                                      *

return

10 format(1X,A)
100 format(I5,1X,I3,1X,I4,1X,F12.6)
101 format(I5,1X,I3,1X,I4,1X,F12.6,1X,'[',I4,1X,F12.6,']')

end subroutine Print_Isotopes
