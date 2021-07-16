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
! Copyright (C) 1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************

subroutine DETSTR2(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,ISGN,IWORK,IPRINT)
! AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
! MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
! PURPOSE:
!
! A DETERMINANT, IDET, IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
! POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
! INDICATES BETA ORBITAL.
! FIND CORRESPONDING ALPHA STRING AND BETA STRING

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NEL, NAEL, NBEL, IDET(NEL), IPRINT
integer(kind=iwp), intent(out) :: IASTR(NAEL), IBSTR(NBEL), ISGN, IWORK(NEL)
integer(kind=iwp) :: IBEL, IEXPO

! FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
! THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING

call ORDSTR(IDET,IWORK,NEL,ISGN,IPRINT)

! ALPHA STRING IS LAST NAEL ORBITALS
call ICOPY(NAEL,IWORK(NBEL+1),1,IASTR,1)

! BETA STRING MUST BE COMPLETELY TURNED AROUND
do IBEL=1,NBEL
  IBSTR(IBEL) = -IWORK(NBEL+1-IBEL)
end do
! SIGN CHANGE FOR SWITCH OF BETA ORBITALS
IEXPO = (NBEL*NBEL+NBEL)/2
ISGN = ISGN*(-1)**IEXPO

if (IPRINT > 30) then
  write(u6,*) ' INPUT DETERMINANT'
  call IWRTMA(IDET,1,NEL,1,NEL)
  write(u6,*) ' CORRESPONDING ALPHA STRING'
  call IWRTMA(IASTR,1,NAEL,1,NAEL)
  write(u6,*) ' CORRESPONDING BETA STRING'
  call IWRTMA(IBSTR,1,NBEL,1,NBEL)
  write(u6,*) ' ISGN FOR SWITCH ',ISGN
end if

return

end subroutine DETSTR2
