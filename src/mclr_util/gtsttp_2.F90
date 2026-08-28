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

subroutine GTSTTP_2(ICLS,IEL1,IEL3,ITYPE)
! Relation between number of electrons in RAS1, RAS3 and string type
!
! GET IEL1,IEL3 : Number of electrons of class ICLS of type ITYPE
!
! Jeppe Olsen, Another lonely night in Lund

use Str_Info, only: ITYP_Dummy, Str
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ICLS, ITYPE
integer(kind=iwp), intent(out) :: IEL1, IEL3

call GTSTTPS_2(IEL1,IEL3,Str(min(ICLS,ITYP_Dummy))%EL1,Str(min(ICLS,ITYP_Dummy))%EL3,ITYPE)

end subroutine GTSTTP_2
