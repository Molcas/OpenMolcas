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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine RWDTG(Num,DMat,lth,Option,DT,iDisk,MaxNum)
!***********************************************************************
!                                                                      *
! purpose: Read / write density matrix, two-electron hamiltonian       *
!          or gradient to the file if they are on disk / can't be      *
!          stored in core                                              *
!                                                                      *
! input:                                                               *
!   Num       - 'map number'                                           *
!   DMat(lth) - matrix to be written or read (Dens or TwoHam)          *
!   Option    - R- read, W- write                                      *
!   DT        - 'DENS  ', 'dVxcdR', 'TWOHAM', 'GRAD  '                 *
!   iDisk     - position in the file                                   *
!   MaxNum    - maximal number of density matrices that can be         *
!               stored on the disk                                     *
!                                                                      *
! output:                                                              *
!   DMat      - matrix read from the disk (if the option was read)     *
!                                                                      *
!***********************************************************************

use Files, only: LuDSt, LuGrd, LuOSt, LuTSt

implicit none
integer Num, lth, MaxNum
real*8 DMat(lth)
integer iDisk(MaxNum)
character Option, DT*6
integer jDisk, LU
#include "SysDef.fh"

! Check Num; Subroutine is called with Num = - MapDns(i)

if (Num <= 0) then
  write(6,*) 'RWDTG: Num <= 0'
  write(6,*) 'Num=',Num
  write(6,*) 'Wrong density number supplied.'
  call Abend()
end if
if (Num > MaxNum) then
  write(6,*) 'RWDTG: Num > MaxNum'
  write(6,*) 'Num,MaxNum=',Num,MaxNum
  write(6,*) 'Wrong density number supplied.'
  call Abend()
end if

if ((DT /= 'DENS') .and. (DT /= 'TWOHAM') .and. (DT /= 'GRAD') .and. (DT /= 'dVxcdR')) then
  write(6,*) 'RWDTG: invalid value of DT'
  write(6,*) '->DT<-=->',DT,'<-'
  write(6,*) 'Valid values: "DENS  "'
  write(6,*) '              "dVxcdR"'
  write(6,*) '              "TWOHAM"'
  write(6,*) '              "GRAD  "'
  call Abend()
end if

if ((Option /= 'W') .and. (Option /= 'R')) then
  write(6,*) 'RWDTG: invalid Option'
  write(6,*) '->Option<-=->',Option,'<-'
  write(6,*) 'Valid Options: R'
  write(6,*) '               W'
end if

if (DT == 'DENS  ') then
  LU = LuDSt
else if (DT == 'TWOHAM') then
  LU = LuTSt
else if (DT == 'GRAD  ') then
  LU = LuGrd
else
  LU = LuOSt
end if

if (Option == 'W') then

  ! Write density matrix to DNS

  if (Num == 1) iDisk(Num) = 0
  jDisk = iDisk(Num)
  if (jDisk == -1) then
    write(6,*) 'RWDTG: jDisk == -1'
    write(6,*) 'Num,MaxNum=',Num,MaxNum
    write(6,*) 'The preceeding block was not written.'
    call Abend()
  end if
  call dDaFile(LU,1,DMat,lth,jDisk)
  if (Num+1 <= MaxNum) iDisk(Num+1) = jDisk

else if (Option == 'R') then

  ! Read density matrix from DNS

  jDisk = iDisk(Num)
  call dDaFile(LU,2,DMat,lth,jDisk)
end if

end subroutine RWDTG
