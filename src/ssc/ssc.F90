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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine ssc(iRC)
use definitions, only: iwp
use spool, only: SpoolInp, Close_LuSpool
Implicit None
integer(kind=iwp), Intent(out) :: iRC

integer(kind=iwp) :: LuSpool, nDiff
logical(kind=iwp) :: DoRys
integer(kind=iwp), external :: IsFreeUnit

LuSpool = 37
LuSpool=IsFreeUnit(LuSpool)
call SpoolInp(LuSpool)

nDiff = 2
DoRys = .True.
call IniSew(DoRys,nDiff)

call Close_LuSpool(LuSpool)

Call Drv2El_BP()

Call ClsSew()
iRC=0

end subroutine ssc
