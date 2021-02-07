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

implicit real*8(a-h,o-z)
#include "MD.fh"
#include "WrkSpc.fh"
!character*180 Lines(10)

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
!write(6,*)
!write(6,'(a,i8,a)')   ' The Dynamix program has ',mxMem,' double precision words memory available.'
!write(6,'(a,i8,a)')   '                       ',mxMem*8/(1024**2),' MB'
!write(6,*)
!#endif

! Set the default values

call Get_nAtoms_Full(natom)
THERMO = 0
TEMP = 298.15
VELO = 0
POUT = 0
PIN = natom*3
iPrint = 2
DT = 1.0d1
RESTART = 0.0d0
lH5Restart = .false.

!Etot0 = 0.0D0

return

end subroutine Init_Dynamix
