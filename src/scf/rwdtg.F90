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

use SCFFiles, only: LuDSt, LuGrd, LuOSt, LuTSt
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Num, lth, MaxNum
real(kind=wp), intent(inout) :: DMat(lth)
character, intent(in) :: Option
integer(kind=iwp), intent(inout) :: iDisk(MaxNum)
character(len=6) :: DT
integer(kind=iwp) :: jDisk, LU

! Check Num; Subroutine is called with Num = - MapDns(i)

if (Num <= 0) then
  write(u6,*) 'RWDTG: Num <= 0'
  write(u6,*) 'Num=',Num
  write(u6,*) 'Wrong density number supplied.'
  call Abend()
end if
if (Num > MaxNum) then
  write(u6,*) 'RWDTG: Num > MaxNum'
  write(u6,*) 'Num,MaxNum=',Num,MaxNum
  write(u6,*) 'Wrong density number supplied.'
  call Abend()
end if

select case (DT)
  case ('DENS')
    LU = LuDSt
  case ('TWOHAM')
    LU = LuTSt
  case ('GRAD')
    LU = LuGrd
  case ('dVxcdR')
    LU = LuOSt
  case default
    write(u6,*) 'RWDTG: invalid value of DT'
    write(u6,*) '->DT<-=->',DT,'<-'
    write(u6,*) 'Valid values: "DENS  "'
    write(u6,*) '              "dVxcdR"'
    write(u6,*) '              "TWOHAM"'
    write(u6,*) '              "GRAD  "'
    call Abend()
end select

select case (Option)

  case ('W')

    ! Write density matrix to DNS

    if (Num == 1) iDisk(Num) = 0
    jDisk = iDisk(Num)
    if (jDisk == -1) then
      write(u6,*) 'RWDTG: jDisk == -1'
      write(u6,*) 'Num,MaxNum=',Num,MaxNum
      write(u6,*) 'The preceeding block was not written.'
      call Abend()
    end if
    call dDaFile(LU,1,DMat,lth,jDisk)
    if (Num+1 <= MaxNum) iDisk(Num+1) = jDisk

  case ('R')

    ! Read density matrix from DNS

    jDisk = iDisk(Num)
    call dDaFile(LU,2,DMat,lth,jDisk)

  case default

    write(u6,*) 'RWDTG: invalid Option'
    write(u6,*) '->Option<-=->',Option,'<-'
    write(u6,*) 'Valid Options: R'
    write(u6,*) '               W'

end select

end subroutine RWDTG
