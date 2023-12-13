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

do
  key = Get_Ln(LuSpool)
  Line = adjustl(key)
  if (Line(1:1) == '*') cycle
  if (Line == ' ') cycle
  call UpCase(Line)
  select case (Line(1:4))
    case ('NOEX')
      DoExpbas = .false.
    case ('DESY')
      DoDesy = .true.
    case ('FILE')
      Line = Get_Ln(LuSpool)
      call FileOrb(Line,EB_FileOrb)
    case ('END ')
      exit
    case default
      write(u6,*) 'Unidentified key word  : '
      call FindErrorLine()
      call Quit_OnUserError()
  end select
end do

! END of Input

end subroutine Readinp_expbas
