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

!#define _DEBUGPRINT_
subroutine ANNSTR(STRING,NSTINI,NSTINO,NEL,NORB,Z,NEWORD,LROW,LSGSTR,ISGSTI,ISGSTO,TI,TTO,I1TYP)
! A set of strings containing NEL electrons are given
! set up all possible ways of annihilating an electron from
! this set of string
!
!========
! Input :
!========
! STRING : Input strings containing NEL electrons
! NSTINI : Number of NEL  strings
! NSTINO : Number of NEL-1 strings
! NEL    : Number of electrons in input strings
! NORB   : Number of orbitals
! Z      : Lexical ordering matrix for output strings containing
!          NEL - 1 electrons
! NEWORD : Reordering array for N-1 strings
! LROW   : Number of rows in output vector
!          = NEL : Information is written in truncated form
!                  row  corresponds to place of electron
!          = NORB: Information is written in expanded form,
!                  row corresponds to full orbital number
! LSGSTI : /= 0 => use ISGSTI,ISGSTO to allow for sign of string
! ISGSTI : Sign array for NEL   strings
! ISGSTO : Sign array for NEL-1 strings
!
!=========
! Output :
!=========
!TI      : Array giving minus orbital annihilated
!          TI(I,ISTRIN) < 0 : Orbital TI(I,ISTRIN) can be
!          annihilated from string ISTRIN
!TTO     : Resulting NEL - 1 strings
!          if the resulting string has carries a negative sign
!          then the string number is shifted with a NSTINO
!          TTO(I,ISTRIN) = <= 0 indicates that orbital I cannot be
!                             added to istrin
!          TTO(I,ISTRIN) = LSTRIN /= 0 indicates that I added to
!                          ISTRIN  gives LSTRIN. If LSTRIN is
!                          nehative
!                          A+(I)!ISTRIN>=-!-LSTRIN>.

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSTINI, NEL, STRING(NEL,NSTINI), NSTINO, NORB, Z(NORB,NEL-1), NEWORD(NSTINO), LROW, LSGSTR, &
                                 ISGSTI(NSTINI), ISGSTO(NSTINO), I1TYP
integer(kind=iwp), intent(out) :: TI(LROW,NSTINI), TTO(LROW,NSTINI)
integer(kind=iwp) :: IEL, IEXPN, IIISGN, IPLACE, ISTRIN, ITYPE, JSTRIN, STRIN2(NEL-1)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I, NPR
integer(kind=iwp), parameter :: MAXPR = 60
#endif
integer(kind=iwp), external :: IOCTP2_MCLR, ISTRNM

#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' ANNSTR speaking'
write(u6,*) ' ==============='
#endif
! Expanded or truncated form
if ((LROW == NEL) .and. (NEL /= NORB)) then
  IEXPN = 0
else
  IEXPN = 1
end if
! Loop over input strings
do ISTRIN=1,NSTINI
  ! loop over electrons to be removed
  do IEL=1,NEL
    if (IEXPN == 0) then
      IPLACE = IEL
    else
      IPLACE = STRING(IEL,ISTRIN)
    end if
    STRIN2(1:IEL-1) = STRING(1:IEL-1,ISTRIN)
    STRIN2(IEL:NEL-1) = STRING(IEL+1:NEL,ISTRIN)
    ! Is new string allowed ?
    ITYPE = IOCTP2_MCLR(STRIN2,NEL-1,I1TYP)
    if (ITYPE /= 0) then
      !        ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
      JSTRIN = ISTRNM(STRIN2,NORB,NEL-1,Z,NEWORD,1)
      TTO(IPLACE,ISTRIN) = JSTRIN
      TI(IPLACE,ISTRIN) = -STRING(IEL,ISTRIN)
      IIISGN = (-1)**(IEL-1)
      if (LSGSTR > 0) IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
      if (IIISGN == -1) TTO(IPLACE,ISTRIN) = -TTO(IPLACE,ISTRIN)
    end if
  end do
end do

#ifdef _DEBUGPRINT_
NPR = min(NSTINI,MAXPR)
write(u6,*) ' Output from ANNSTR :'
write(u6,*) '==================='
if (IEXPN == 0) then
  write(u6,*) ' Strings with an electron removed'
else
  write(u6,*) ' Combined N+1/N-1 string array'
end if
do ISTRIN=1,NPR
  write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' New strings.. ',(TTO(I,ISTRIN),I=1,LROW)
end do

write(u6,*) ' orbitals removed'
do ISTRIN=1,NPR
  write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' orbitals annihilated.. ',(TI(I,ISTRIN),I=1,LROW)
end do
#endif

return

end subroutine ANNSTR
