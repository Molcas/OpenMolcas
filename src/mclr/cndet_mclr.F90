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

!#define _DEBUGPRINT_
subroutine CNDET_MCLR(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,IDET)
! A configuration ICONF in compressed form and a set of
! prototype determinants, IPDET, are given.
!
! Construct the corresponding determinants in contracted form.
!
! JEPPE OLSEN, NOVEMBER 1988

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB, ICONF(NORB), NDET, NOP, IPDET(NOP,NDET), NEL, NCL
integer(kind=iwp), intent(out) :: IDET(NEL,NDET)
integer(kind=iwp) :: IADD, IBASE, ICL, IOP, JDET

! POSITIVE NUMBER: ALPHA ORBITAL
! NEGATIVE NUMBER: BETA  ORBITAL

#ifdef _DEBUGPRINT_
if (NCL /= 0) then
  write(u6,*) ' DOUBLE OCCUPIED ORBITALS'
  call IWRTMA(ICONF,1,NCL,1,NCL)
end if
if (NOP /= 0) then
  write(u6,*) ' OPEN ORBITALS'
  call IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
end if
#endif

!1  DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST

do ICL=1,NCL
  IBASE = 2*(ICL-1)
  IDET(IBASE+1,:) = ICONF(ICL)
  IDET(IBASE+2,:) = -ICONF(ICL)
end do

!2  SINGLY OCCUPIED ORBITALS

IADD = 2*NCL
do JDET=1,NDET
  do IOP=1,NOP
    if (IPDET(IOP,JDET) == 1) IDET(IADD+IOP,JDET) = ICONF(NCL+IOP)
    if (IPDET(IOP,JDET) == 0) IDET(IADD+IOP,JDET) = -ICONF(NCL+IOP)
  end do
end do

!3  OUTPUT

#ifdef _DEBUGPRINT_
write(u6,*) ' CONFIGURATION FROM DETCON'
call IWRTMA(ICONF,1,NORB,1,NORB)
if (NOP*NDET > 0) then
  write(u6,*) ' PROTO TYPE DETERMINANTS'
  call IWRTMA(IPDET,NOP,NDET,NOP,NDET)
end if
if (NEL*NDET > 0) then
  write(u6,*) ' CORRESPONDING DETERMINANTS'
  call IWRTMA(IDET,NEL,NDET,NEL,NDET)
end if
#endif

end subroutine CNDET_MCLR
