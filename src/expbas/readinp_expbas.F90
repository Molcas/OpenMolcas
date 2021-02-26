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
! Copyright (C) Giovanni Li Manni                                      *
!***********************************************************************

subroutine Readinp_expbas()
! Author: G. Li Manni (University of Geneva)

use info_expbas_mod, only: DoExpbas, DoDesy, EB_FileOrb
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: LuSpool
character(len=180) :: Line, key
integer(kind=iwp), external :: isFreeUnit
character(len=180), external :: Get_Ln

! Initial values

DoExpbas = .true.
DoDesy = .false.
EB_FileOrb = ' '

LuSpool = 18
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'EXPBAS')

999 continue
!read(LuSpool,'(A)',End=9940) Line
key = Get_Ln(LuSpool)
call LeftAd(key)
Line = key
if (Line(1:1) == '*') goto 999
if (Line == ' ') goto 999
call UpCase(Line)
if (Line(1:4) == 'NOEX') goto 1000
if (Line(1:4) == 'DESY') goto 2000
if (Line(1:4) == 'FILE') goto 3000
if (Line(1:4) == 'END ') Go To 99999
write(u6,*) 'Unidentified key word  : '
call FindErrorLine()
call Quit_OnUserError()

!========= NOEX =============
1000 continue
DoExpbas = .false.
Go To 999

!========= DESY =============
2000 continue
DoDesy = .true.
Go To 999

!========= FILE =============
3000 continue
Line = Get_Ln(LuSpool)
call FileOrb(Line,EB_FileOrb)
Go To 999

! END of Input

!9940  Continue
write(u6,*) ' READIN: Premature end of file when reading selected'
call Abend()

99999 continue

end subroutine Readinp_expbas
