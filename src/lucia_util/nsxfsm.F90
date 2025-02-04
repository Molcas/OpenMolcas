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

function NSXFSM(NSMOB,MXPOBS,NO1PS,NO2PS,ISXSM,ADSXA,ISYM,IPRNT)
! Number of single excitations of symmetry ISXSM
!
! ISYM = 0 : All symmetry allowed excitations
! ISYM = 1 : Only excitations a+iaj with I >= J
! ISYM =-1 : Only excitations a+iaj with I > J

use Definitions, only: u6

integer ADSXA(MXPOBS,2*MXPOBS)
integer NO1PS(*), NO2PS(*)

MSXFSM = 0
!write(u6,*) ' NSMOB ',NSMOB
do IO1SM=1,NSMOB
  IO2SM = ADSXA(IO1SM,ISXSM)
  !write(u6,*) ' IO1SM,IO2SM',IO1SM,IO2SM
  if ((ISYM == 0) .or. (IO1SM > IO2SM)) then
    MSXFSM = MSXFSM+NO1PS(IO1SM)*NO2PS(IO2SM)
  else if ((ISYM == 1) .and. (IO1SM == IO2SM)) then
    MSXFSM = MSXFSM+NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
  else if ((ISYM == -1) .and. (IO1SM == IO2SM)) then
    MSXFSM = MSXFSM+NO1PS(IO1SM)*(NO1PS(IO1SM)-1)/2
  end if
end do

NSXFSM = MSXFSM

NTEST = 0
NTEST = max(NTEST,IPRNT)

if (NTEST /= 0) write(u6,*) ' Number of single excitations of symmetry ',ISXSM,',',NSXFSM

end function NSXFSM
