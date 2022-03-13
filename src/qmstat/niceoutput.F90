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

subroutine NiceOutPut(EelP,Gam,Gamma,BetaBol)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "real.fh"
#include "warnings.h"
#include "constants.fh"
parameter(Conver1=1.0d10*CONST_BOHR_RADIUS_IN_SI_)
parameter(Conver2=2.0d0*Pi/360.0d0)
character*3 EelP
character*40 Word1, Word2, Word3
logical Eq, Pr, It, Cl, Qu
external Len_TrimAO

! Enter.

! Make sure that the string is of correct length.

if (Len_TrimAO(EelP) /= 3) then
  write(6,*) 'Illegal call to NiceOutPut'
  call Quit(_RC_INTERNAL_ERROR_)
end if

! Check what type of output that is requested.

Eq = .false.
Pr = .false.
It = .false.
Cl = .false.
Qu = .false.
if (index(EelP,'E') /= 0) Eq = .true.
if (index(EelP,'P') /= 0) Pr = .true.
if (index(EelP,'I') /= 0) It = .true.
if (index(EelP,'C') /= 0) Cl = .true.
if (index(EelP,'Q') /= 0) Qu = .true.

! Start printing!

! With aid of concatenation, we here construct a header.

write(6,*)
write(6,*)
write(6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(6,*) '  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  '
write(6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(6,*)
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
iM1 = Len_TrimAO(Word1)
iM2 = Len_TrimAO(Word2)
iM3 = Len_TrimAO(Word3)
write(6,*) Word1(1:iM1)//Word2(1:iM2)//Word3(1:iM3)

! Now dump a lot of information.

if (It) then
  write(6,*)
  write(6,*)
  write(6,11) '*  Parameters of the calculation  *'
  write(6,*)
  write(6,12) '--Macroscopic quantities'
  write(6,12) '  Temperature(K)      Pressure(Atm.)   Permitivity'
  write(6,13) Temp,Pres,Diel
  write(6,12) '--Maximal MC-Step parameters'
  write(6,12) '  Translation(Ang.)   Rotation(deg.)   Cavity Radius(Ang.)'
  write(6,13) delX*Conver1,delFi/Conver2,delR*Conver1
  write(6,12) '--Configuration data'
  write(6,12) '  Initial conf.       Writing conf.    MC-Steps'
  if (iNrIn >= 0) write(6,14) iNrIn,iNrUt,nMicro*nMacro
  if (iNrIn < 0) write(6,15) '   Random/Input',iNrUt,nMicro*nMacro
  if (QmType(1:4) /= 'RASS') then
    write(6,12) '--Hartree-Fock simulation data'
    write(6,12) '  Total Occupation    Number of Orbitals'
    write(6,16) iOcc1,iOrb(1)
  else
    write(6,12) '--Rassi state simulation data'
    write(6,12) '  State interacting with solvent'
    if (.not. lCiSelect) then
      write(6,17) nEqState
    else
      write(6,12) '   CI-select overlap option used'
    end if
    write(6,12) '  State threshold     Density threshold'
    if (MoAveRed .and. ContrStateB) then
      write(6,18) ThrsCont,ThrsRedOcc
    else if (MoAveRed .and. (.not. ContrStateB)) then
      write(6,19) '   N/A',ThrsRedOcc
    else if ((.not. MoAveRed) .and. ContrStateB) then
      write(6,20) ThrsCont,'N/A'
    else
      write(6,12) '   N/A                 N/A'
    end if
    if (nLvlShift /= 0) then
      write(6,12) '  Level shift applied'
    end if
  end if
end if
write(6,*)
write(6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(6,*) '  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  '
write(6,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(6,*)
if (iT) then
  write(6,*)
  write(6,*) 'Simulation progress.'
  write(6,*)
end if

! Some formats

11 format('               ',A)
12 format('    ',A)
13 format('    ',3(F10.4,'        '))
14 format('    ',3(I8,'        '))
15 format('    ',A,2(I8,'         '))
16 format('    ',2(I5,'                '))
17 format('    ',I5)
18 format('     ',2(E11.4,'         '),' ')
19 format('    ',A,'               ',E11.4)
20 format('     ',E11.4,'           ',A)

! Tschuss

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(Gam)
  call Unused_real(Gamma)
  call Unused_real(BetaBol)
end if

end subroutine NiceOutPut
