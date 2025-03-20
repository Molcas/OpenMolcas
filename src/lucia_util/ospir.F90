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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine OSPIR(NOSPIR,IOSPIR,MXPIRR,MXPOBS)
! Number and symmetries of orbitals corresponding to a given shell
!
! =====
! Input
! =====
!
! PNTGRP  : type of pointgroup
!       = 1 => D2h or a subgroup of D2H
!       = 2 => C inf v
!       = 3 => D inf h
!       = 4 => O 3
! MXPIRR : Largest allowed number of shell irreps
! MXPOBS : Largest allowed number of orbital symmetries
!
! ======
! Output
! ======
!
! NOSPIR : Number of orbital symmetries per irrep
! IOSPIR : Orbital symmetries corresponding to a given irrep
!
! Jeppe Olsen, Winter of 1991

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use lucia_data, only: NIRREP
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MXPIRR, MXPOBS
integer(kind=iwp), intent(out) :: NOSPIR(MXPIRR), IOSPIR(MXPOBS,MXPIRR)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IRREP
#endif

!if (PNTGRP == 1) then
!=====
! D2h
!=====
NOSPIR(1:8) = 1
IOSPIR(1,1:8) = [1,2,3,4,5,6,7,8]
!else
!  write(u6,*) ' Sorry  PNTGRP out of range, PNTGRP = ',PNTGRP
!  write(u6,*) ' OSPIR fatally wounded'
!  !stop 5
!  call SYSABENDMSG('lucia_util/ospir','Internal error','')
!end if

#ifdef _DEBUGPRINT_
write(u6,*) ' OSPIR speaking'
write(u6,*) ' ================'
write(u6,*) ' Number of orbitals per irrep'
call IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
write(u6,*) ' Orbital symmetries per irrep'
do IRREP=1,NIRREP
  call IWRTMA(IOSPIR(:,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
end do
#endif

end subroutine OSPIR
