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

subroutine CNDET(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,IDET,IPRINT)
!     AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
!     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
!     PURPOSE:
!
!     A CONFIGURATION ICONF IN COMPRESSED FORM AND A SET OF
!     PROTOTYPE DETERMINANTS ,IPDET, ARE GIVEN .
!     CONSTRUCT THE CORRESPONDING DETERMINANTS IN CONTRACTED  FORM .
!
!     SUBROUTINE CALLS:
!
!     IWRTMA

implicit real*8(A-H,O-Z)
dimension ICONF(*)
dimension IPDET(*)
dimension IDET(NEL,*)

! POSITIVE NUMBER  : ALPHA ORBITAL
! NEGATIVE NUMBER  : BETA  ORBITAL
!
!IPRINT = 40
!write(6,*) 'iprint == 40 in CNDET, nclosed, nopen',iprint, ncl,nop
!write(6,*) 'iconf, nclosed, nopen',iconf(1:1+ncl-1)
!if (nop /= 0) then
!  write(6,*) 'opened orbitals',iconf(1+ncl:ncl+nop)
!end if

if (IPRINT == 40) then
  if (NCL /= 0) then
    write(6,*) ' DOUBLE OCCUPIED ORBITALS '
    call IWRTMA(ICONF,1,NCL,1,NCL)
  end if
  if (NOP /= 0) then
    write(6,*) ' OPEN ORBITALS '
    call IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
  end if
end if

! 1 DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST

do ICL=1,NCL
  IBASE = 2*(ICL-1)
  do JDET=1,NDET
    IDET(IBASE+1,JDET) = ICONF(ICL)
    IDET(IBASE+2,JDET) = -ICONF(ICL)
  end do
end do

! 2 SINGLY OCCUPIED ORBITALS

IADD = 2*NCL
do JDET=1,NDET
  do IOP=1,NOP
    IADR = (JDET-1)*NOP+IOP
    if (IPDET(IADR) == 1) IDET(IADD+IOP,JDET) = ICONF(NCL+IOP)

    if (IPDET(IADR) == 0) IDET(IADD+IOP,JDET) = -ICONF(NCL+IOP)

  end do
end do

if (IPRINT == 40) then
  write(6,*) ' CONFIGURATION FROM DETCON '
  call IWRTMA(ICONF,1,NORB,1,NORB)
  write(6,*) ' PROTO TYPE DETERMINANTS '
  if (NOP*NDET > 0) call IWRTMA(IPDET,NOP,NDET,NOP,NDET)

  if (NEL*NDET > 0) write(6,*) ' CORRESPONDING DETERMINANTS '
  call IWRTMA(IDET,NEL,NDET,NEL,NDET)

end if

IPRINT = 0

return

end subroutine CNDET
