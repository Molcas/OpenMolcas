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
subroutine SHTOOB(NSHPIR,NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,NOBPS,NOB)
! Number of shells per irrep => Number of orbitals per symmetry
!
! =====
! Input
! =====
!  NSHPIR : Number of shells per irrep
!  NIRREP : Number of irreps
!  MXPOBS : Largest allowed number of orbitals symmetries
!  NSMOB  : Number of orbital symmetries
!  NOSPIR : Number of orbital symmetries per irrep
!  IOSPIR : Orbital symmetries per irrep
!
! ======
! Output
! ======
!  NOBPS  : Number of orbitals per symmetry
!  NOB    : Number of orbitals
!
! Jeppe Olsen, Winter of 1991

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NIRREP, NSHPIR(NIRREP), MXPOBS, NSMOB, NOSPIR(NIRREP), IOSPIR(MXPOBS,NIRREP)
integer(kind=iwp), intent(out) :: NOBPS(NSMOB), NOB
integer(kind=iwp) :: IISM, IRREP, ISM

NOBPS(1:NSMOB) = 0
NOB = 0
do IRREP=1,NIRREP
  do ISM=1,NOSPIR(IRREP)
    IISM = IOSPIR(ISM,IRREP)
    NOBPS(IISM) = NOBPS(IISM)+NSHPIR(IRREP)
    NOB = NOB+NSHPIR(IRREP)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' SHTOOB Speaking'
write(u6,*) ' ==============='
write(u6,*) ' Number of orbitals obtained ',NOB
write(u6,*) ' Number of orbitals per symmetry'
call IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
#endif

end subroutine SHTOOB
