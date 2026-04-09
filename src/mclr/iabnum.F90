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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

function IABNUM(IASTR,IBSTR,IAGRP,IBGRP,IGENSG,ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN)
! Encapsulation routine for IABNUS

use Str_info, only: nElec, NoCTyp, STR
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IABNUM
integer(kind=iwp), intent(in) :: IASTR(*), IBSTR(*), IAGRP, IBGRP, IGENSG, ISGNA(*), ISGNB(*), &
                                 IOOS(NOCTYP(IAGRP),NOCTYP(IBGRP),*), NORB
integer(kind=iwp), intent(out) :: ISGNAB, IPSFAC
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp), external :: IABNUS

IABNUM = IABNUS(IASTR,NELEC(IAGRP),Str(IAGRP)%STREO,Str(IAGRP)%STCL,Str(IAGRP)%STSM,NOCTYP(IAGRP),Str(IAGRP)%Z,Str(IAGRP)%ISTSO, &
                Str(IAGRP)%NSTSO,IBSTR,NELEC(IBGRP),Str(IBGRP)%STREO,Str(IBGRP)%STCL,Str(IBGRP)%STSM,NOCTYP(IBGRP),Str(IBGRP)%Z, &
                Str(IBGRP)%ISTSO,Str(IBGRP)%NSTSO,IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,PSSIGN,IPSFAC)

end function IABNUM
