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
! prototype determinants,IPDET, are given.
!
! Construct the corresponding determinants in contracted  form.
!
! JEPPE OLSEN, NOVEMBER 1988

implicit none
integer NEL
integer ICONF(*)
!integer IPDET(NOP,NDET)
integer IPDET(*)
integer NORB, NOP, NCL
integer IDET(NEL,*)
! local variables
integer ICL, IBASE, JDET, NDET, IADD, IOP, IADR

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(NORB)
#endif

! POSITIVE NUMBER: ALPHA ORBITAL
! NEGATIVE NUMBER: BETA  ORBITAL

#ifdef _DEBUGPRINT_
if (NCL /= 0) then
  write(6,*) ' DOUBLE OCCUPIED ORBITALS'
  call IWRTMA(ICONF,1,NCL,1,NCL)
end if
if (NOP /= 0) then
  write(6,*) ' OPEN ORBITALS'
  call IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
end if
#endif

!1  DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST

do ICL=1,NCL
  IBASE = 2*(ICL-1)
  do JDET=1,NDET
    IDET(IBASE+1,JDET) = ICONF(ICL)
    IDET(IBASE+2,JDET) = -ICONF(ICL)
  end do
end do

!2  SINGLY OCCUPIED ORBITALS

IADD = 2*NCL
do JDET=1,NDET
  do IOP=1,NOP
    IADR = (JDET-1)*NOP+IOP
    if (IPDET(IADR) == 1) IDET(IADD+IOP,JDET) = ICONF(NCL+IOP)
    if (IPDET(IADR) == 0) IDET(IADD+IOP,JDET) = -ICONF(NCL+IOP)
  end do
end do

!3  OUTPUT

#ifdef _DEBUGPRINT_
write(6,*) ' CONFIGURATION FROM DETCON'
call IWRTMA(ICONF,1,NORB,1,NORB)
write(6,*) ' PROTO TYPE DETERMINANTS'
if (NOP*NDET > 0) call IWRTMA(IPDET,NOP,NDET,NOP,NDET)
if (NEL*NDET > 0) write(6,*) ' CORRESPONDING DETERMINANTS'
call IWRTMA(IDET,NEL,NDET,NEL,NDET)
#endif

end subroutine CNDET_MCLR
