!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ANNSTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NACOB,IEC,LDIM,IPRNT)
! A group of strings containing NEL electrons is given
! Set up all possible ways of removing an electron
!
!========
! Input :
!========
! STRING : Input strings containing NEL electrons
! NSTINI : Number of input  strings
! NSTINO : Number of output strings
! NEL    : Number of electrons in input strings
! NORB   : Number of orbitals
! IORBOF : Number of first orbital
! Z      : Lexical ordering matrix for output strings containing
!          NEL - 1 electrons
! NEWORD : Reordering array for N-1 strings
! LSGSTR : /= 0 => Include sign arrays ISGSTI,ISGSTO of strings
! ISGSTI : Sign array for NEL   strings
! ISGSTO : Sign array for NEL-1 strings
! IEC    : = 1 Extended map, dimension equals number of orbs
! IEC    : = 2 Compact  map, dimension equals number of elecs
! LDIM   : Row dimension ( see IEC)
!
!=========
! Output :
!=========
!
!TI      : TI(I,ISTRIN) > 0 indicates that orbital I can be added
!          to string ISTRIN .
!TTO     : Resulting NEL + 1 strings
!          if the string have a negative sign
!          then the phase equals - 1

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NSTINI, NEL, STRING(NEL,NSTINI), NSTINO, NORB, IORBOF, Z(NORB,NEL+1), NEWORD(NSTINO), LSGSTR, ISGSTI(NSTINI), &
                     ISGSTO(NSTINO), LDIM, TI(LDIM,NSTINI), TTO(LDIM,NSTINI), NACOB, IEC, IPRNT
integer(kind=iwp) :: I, IEL, IIISGN, IORB, IORBABS, IROW, ISTRIN, JEL, JSTRIN, MAXPR, NPR, NTEST, NTEST0, STRIN2(500)
integer(kind=iwp), external :: ISTRNM

NTEST0 = 1
NTEST = max(IPRNT,NTEST0)
if (NTEST >= 20) then
  write(u6,*) ' ==============='
  write(u6,*) ' ANNSTR speaking'
  write(u6,*) ' ==============='
  write(u6,*)
  write(u6,*) ' Number of input electrons ',NEL
end if

do ISTRIN=1,NSTINI
  do IEL=1,NEL
    ! String with electron removed
    do JEL=1,IEL-1
      STRIN2(JEL) = STRING(JEL,ISTRIN)
    end do
    do JEL=IEL+1,NEL
      STRIN2(JEL-1) = STRING(JEL,ISTRIN)
    end do
    JSTRIN = ISTRNM(STRIN2,NACOB,NEL-1,Z,NEWORD,1)
    !write(u6,*) ' anni-string and number'
    !call IWRTMA(STRIN2,1,NEL-1,1,NEL-1)
    !write(u6,*) ' JSTRIN = ',JSTRIN

    IORBABS = STRING(IEL,ISTRIN)
    IORB = STRING(IEL,ISTRIN)-IORBOF+1
    if (IEC == 1) then
      IROW = IORB
    else
      IROW = IEL
    end if

    TI(IROW,ISTRIN) = -IORBABS
    !TI(IROW,ISTRIN) = -IORB
    TTO(IROW,ISTRIN) = JSTRIN
    !PAM2009 IIISGN = (-1)**(IEL-1)
    IIISGN = 1-2*mod(IEL+1,2)
    if (LSGSTR /= 0) IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
    if (IIISGN == -1) TTO(IROW,ISTRIN) = -TTO(IROW,ISTRIN)
  end do

end do

if (NTEST >= 20) then
  MAXPR = 60
  NPR = min(NSTINI,MAXPR)
  write(u6,*) 'Output from ANNSTR:'
  write(u6,*) '==================='

  write(u6,*)
  write(u6,*) ' Strings with an electron added or removed'
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' New strings.. ',(TTO(I,ISTRIN),I=1,LDIM)
  end do
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' orbitals added or removed ',(TI(I,ISTRIN),I=1,LDIM)
  end do
end if

end subroutine ANNSTR_GAS
