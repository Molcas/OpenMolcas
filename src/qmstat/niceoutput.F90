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

subroutine NiceOutPut(EelP)

use qmstat_global, only: ContrStateB, DelFi, DelR, DelX, Diel, iNrIn, iNrUt, iOcc1, iOrb, lCiSelect, MoAveRed, nEqState, &
                         nLvlShift, nMacro, nMicro, Pres, QmType, Temp, ThrsCont, ThrsRedOcc
use Constants, only: Angstrom, deg2rad
use Definitions, only: iwp, u6

implicit none
character(len=3), intent(in) :: EelP
logical(kind=iwp) :: Cl, Eq, It, Pr, Qu
character(len=40) :: Word1, Word2, Word3
#include "warnings.h"

! Enter.

! Make sure that the string is of correct length.

if (len_trim(EelP) /= 3) then
  write(u6,*) 'Illegal call to NiceOutPut'
  call Quit(_RC_INTERNAL_ERROR_)
end if

! Check what type of output that is requested.

Eq = index(EelP,'E') /= 0
Pr = index(EelP,'P') /= 0
It = index(EelP,'I') /= 0
Cl = index(EelP,'C') /= 0
Qu = index(EelP,'Q') /= 0

! Start printing!

! With aid of concatenation, we here construct a header.

write(u6,*)
write(u6,*)
write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,*) '  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  '
write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,*)
if (It) then
  write(Word1,*) 'QMStat simulation commencing: '
else
  write(Word1,*) 'SampFile analysis commencing'
end if
if (Cl) then
  write(Word2,*) 'All Classical '
else if (Qu) then
  write(Word2,*) 'Combined Quantum-Classical '
else
  write(Word2,*) ' '
end if
if (Eq) then
  write(Word3,*) 'Equilibration'
else if (Pr) then
  write(Word3,*) 'Production'
else
  write(Word3,*) ' '
end if
write(u6,*) trim(Word1)//trim(Word2)//trim(Word3)

! Now dump a lot of information.

if (It) then
  write(u6,*)
  write(u6,*)
  write(u6,11) '*  Parameters of the calculation  *'
  write(u6,*)
  write(u6,12) '--Macroscopic quantities'
  write(u6,12) '  Temperature(K)      Pressure(Atm.)   Permitivity'
  write(u6,13) Temp,Pres,Diel
  write(u6,12) '--Maximal MC-Step parameters'
  write(u6,12) '  Translation(Ang.)   Rotation(deg.)   Cavity Radius(Ang.)'
  write(u6,13) delX*Angstrom,delFi/deg2rad,delR*Angstrom
  write(u6,12) '--Configuration data'
  write(u6,12) '  Initial conf.       Writing conf.    MC-Steps'
  if (iNrIn >= 0) write(u6,14) iNrIn,iNrUt,nMicro*nMacro
  if (iNrIn < 0) write(u6,15) '   Random/Input',iNrUt,nMicro*nMacro
  if (QmType(1:4) /= 'RASS') then
    write(u6,12) '--Hartree-Fock simulation data'
    write(u6,12) '  Total Occupation    Number of Orbitals'
    write(u6,16) iOcc1,iOrb(1)
  else
    write(u6,12) '--Rassi state simulation data'
    write(u6,12) '  State interacting with solvent'
    if (.not. lCiSelect) then
      write(u6,17) nEqState
    else
      write(u6,12) '   CI-select overlap option used'
    end if
    write(u6,12) '  State threshold     Density threshold'
    if (MoAveRed .and. ContrStateB) then
      write(u6,18) ThrsCont,ThrsRedOcc
    else if (MoAveRed .and. (.not. ContrStateB)) then
      write(u6,19) '   N/A',ThrsRedOcc
    else if ((.not. MoAveRed) .and. ContrStateB) then
      write(u6,20) ThrsCont,'N/A'
    else
      write(u6,12) '   N/A                 N/A'
    end if
    if (nLvlShift > 0) write(u6,12) '  Level shift applied'
  end if
end if
write(u6,*)
write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,*) '  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  '
write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,*)
if (iT) then
  write(u6,*)
  write(u6,*) 'Simulation progress.'
  write(u6,*)
end if

! Tschuss

return

11 format('               ',A)
12 format('    ',A)
13 format('    ',3(F10.4,'        '))
14 format('    ',3(I8,'        '))
15 format('    ',A,2(I8,'         '))
16 format('    ',2(I5,'                '))
17 format('    ',I5)
18 format('     ',2(ES11.4,'         '),' ')
19 format('    ',A,'               ',ES11.4)
20 format('     ',ES11.4,'           ',A)

end subroutine NiceOutPut
