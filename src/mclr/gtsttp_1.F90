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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************

subroutine GTSTTP_1(ICLS,IEL1,IEL3,ITYPE)
! Relation between number of electrons in RAS1, RAS3 and string type
!
! Get ITYPE : type of strings belonging to class ICLS
!             with IEL1,IEL3 electrons
!
! Jeppe Olsen, Another lonely night in Lund

use Str_Info, only: ITYP_Dummy, NOCTYP, Str
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ICLS, IEL1, IEL3
integer(kind=iwp), intent(out) :: ITYPE

call GTSTTPS_1(IEL1,IEL3,Str(min(ICLS,ITYP_Dummy))%EL1,Str(min(ICLS,ITYP_Dummy))%EL3,NOCTYP(min(ICLS,ITYP_Dummy)),ITYPE)

end subroutine GTSTTP_1
