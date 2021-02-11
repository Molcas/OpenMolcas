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

subroutine Init_Dynamix()

use Dynamix_globals, only: DT, iPrint, lH5Restart, PIN, POUT, RESTART, TEMP, THERMO, VELO
use Constants, only: Zero, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: natom
!integer(kind=iwp) :: lLine
!character(len=180) :: Lines(10)

!#ifdef _SKIP_
!Lines(1) = _MOLCAS_VERSION_
!#ifdef _DEMO_
!Lines(2) = 'DEMO VERSION '
!#else
!Lines(2) = ' '
!#endif
!Lines(3) = ' '
!Lines(4) = 'D Y N A M I X'
!Lines(5) = ' '
!Lines(6) = 'Written by Igor Schapiro'
!Lines(7) = ' '
!Lines(8) = 'A program for molecular dynamics simulations'
!Lines(9) = ' '
!Lines(10) = 'Code compiled :'//_BUILD_DATE_
!lLine = Len(Lines(1))
!call Banner(Lines,10,lLine)
!
!write(u6,*)
!write(u6,'(a,i8,a)')   ' The Dynamix program has ',mxMem,' double precision words memory available.'
!write(u6,'(a,i8,a)')   '                         ',mxMem*8/(1024**2),' MB'
!write(u6,*)
!#endif

! Set the default values

call Get_nAtoms_Full(natom)
THERMO = 0
TEMP = 298.15_wp
VELO = 0
POUT = 0
PIN = natom*3
iPrint = 2
DT = Ten
RESTART = Zero
lH5Restart = .false.

!Etot0 = Zero

return

end subroutine Init_Dynamix
