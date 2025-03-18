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
! Copyright (C) 1993, Jeppe Olsen                                      *
!***********************************************************************

integer function ISTRN_MCLR(STRING,IGROUP)
! A string belonging to group IGROUP is given.
! find actual number
!
! Jeppe Olsen, September 1993

use Str_Info, only: STR, NELEC
use MCLR_Data, only: NACOB

implicit none
! Specific input
integer STRING(*)
integer IGROUP
integer NEL
integer, external :: ISTRNM

NEL = NELEC(IGROUP)
ISTRN_MCLR = ISTRNM(STRING,NACOB,NEL,Str(IGROUP)%Z,Str(IGROUP)%STREO,1)

end function ISTRN_MCLR
