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
subroutine CRESTR(STRING,NSTINI,NSTINO,NEL,NORB,Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,ISTMPL,ISTMPO,LROW,I1TYP)
! A set of strings containing NEL electrons are given
! set up all possible ways of adding an electron to this set of strings
!
!========
! Input :
!========
!
! STRING : Input strings containing NEL electrons
! NSTINI : Number of input  strings
! NSTINO : Number of output strings
! NEL    : Number of electrons in input strings
! NORB   : Number of orbitals
! Z      : Lexical ordering matrix for output strings containing
!          NEL + 1 electrons
! NEWORD : Reordering array for N+1 strings
! LSGSTR : /= 0 => Include sign arrays ISGSTI,ISGSTO of strings
! ISGSTI : Sign array for NEL   strings
! ISGSTO : Sign array for NEL+1 strings
!
! LROW   : Length of tables for each string, negative number
!          indicates compact form
!=========
! Output :
!=========
!
!TI      : TI(I,ISTRIN) > 0 indicates that orbital I can be added
!          to string ISTRIN.
!TTO     : Resulting NEL + 1 strings
!          if the resulting string has carries a negative sign
!          then the string number is negative
!          TTO(I,ISTRIN) = LSTRIN /= 0 indicates that I added to
!                          ISTRIN  gives LSTRIN. If LSTRIN is
!                          gt. NSTINO then
!                          A+(I)!ISTRIN>=-!LSTRIN-NSTINO>.
! If lrow <= 0 TI and TTO are constructed in sparse form using
! ISTMPO and ISTMPL as pointers
! ISTST

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSTINI, NEL, STRING(NEL,NSTINI), NSTINO, NORB, Z(NORB,NEL+1), NEWORD(NSTINO), LSGSTR, &
                                 ISGSTI(NSTINI), ISGSTO(NSTINO), LROW, I1TYP
integer(kind=iwp), intent(_OUT_) :: TI(*), TTO(*), ISTMPL(*), ISTMPO(*)
integer(kind=iwp) :: IEL, IIISGN, IOFF, IORB, IPLACE, ISTRIN, ITYPE, JSTRIN, LCR, STRIN2(NEL+1)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I, NPR
integer(kind=iwp), parameter :: MAXPR = 60
#endif
integer(kind=iwp), external :: IOCTP2_MCLR, ISTRNM

! LCR NOT DECLARED !!!!!!
!  I SET IT TO ZERO
LCR = 0
#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' CRESTR speaking'
write(u6,*) ' ==============='
#endif

IOFF = 0     ! dummy initialize
IPLACE = -1  ! dummy initialize
do ISTRIN=1,NSTINI
  ! Offset for current creations from this string
  if (ISTRIN == 1) then
    IOFF = 1
  else
    if (LROW > 0) then
      IOFF = (ISTRIN-1)*LROW+1
    else
      IOFF = IOFF+LCR
    end if
  end if
  LCR = 0
  do IORB=1,NORB

    if (NEL == 0) then
      IPLACE = 1

    else if (NEL /= 0) then

      do IEL=1,NEL
        if ((IEL == 1) .and. (STRING(1,ISTRIN) > IORB)) then
          IPLACE = 1
          exit
        else if (((IEL == NEL) .and. (IORB > STRING(IEL,ISTRIN))) .or. &
                 ((IEL < NEL) .and. (IORB > STRING(IEL,ISTRIN)) .and. (IORB < STRING(min(NEL,IEL+1),ISTRIN)))) then
          IPLACE = IEL+1
          exit
        else if (STRING(IEL,ISTRIN) == IORB) then
          IPLACE = 0
          exit
        end if
      end do

    end if

    if (IPLACE == 0) cycle

    ! Generate next string
    STRIN2(1:IPLACE-1) = STRING(1:IPLACE-1,ISTRIN)
    STRIN2(IPLACE) = IORB
    STRIN2(IPLACE+1:NEL+1) = STRING(IPLACE:NEL,ISTRIN)
    ! Is new string allowed?
    ITYPE = IOCTP2_MCLR(STRIN2,NEL+1,I1TYP)
    if (ITYPE /= 0) then
      JSTRIN = ISTRNM(STRIN2,NORB,NEL+1,Z,NEWORD,1)
      ! Save it!
      if (LROW > 0) then
        LCR = IORB
      else
        LCR = LCR+1
      end if

      TTO(IOFF-1+LCR) = JSTRIN
      IIISGN = (-1)**(IPLACE-1)
      if (LSGSTR /= 0) IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
      if (IIISGN == -1) TTO(IOFF-1+LCR) = -TTO(IOFF-1+LCR)
      TI(IOFF-1+LCR) = IORB

      !TTO(IORB,ISTRIN) = JSTRIN
      !IIISGN = (-1)**(IPLACE-1)
      !if (LSGSTR /= 0) IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
      !if (IIISGN == -1) TTO(IORB,ISTRIN) = -TTO(IORB,ISTRIN)
      !TI(IORB,ISTRIN) = IORB
    end if
  end do

  if (LROW < 0) then
    ISTMPO(ISTRIN) = IOFF
    ISTMPL(ISTRIN) = LCR
    !write(u6,*) ' ISTRIN ISTMPO ISTMPL'
    !write(u6,*) ISTRIN,ISTMPO(ISTRIN),ISTMPL(ISTRIN)
  end if
end do

#ifdef _DEBUGPRINT_
NPR = min(NSTINI,MAXPR)
write(u6,*) ' Output from CRESTR :'
write(u6,*) '==================='

if (LROW > 0) then
  write(u6,*) ' Full map'
  write(u6,*)
  write(u6,*) ' Strings with an electron added'
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' New strings.. ',(TTO((ISTRIN-1)*LROW+I),I=1,NORB)
  end do
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' orbitals added or removed ',(TI((ISTRIN-1)*LROW+I),I=1,NORB)
  end do
else
  write(u6,*) ' Compact map'
  write(u6,*)
  write(u6,*) ' Strings with an electron added'
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' New strings.. ',(TTO(ISTMPO(ISTRIN)-1+I),I=1,ISTMPL(ISTRIN))
  end do
  do ISTRIN=1,NPR
    write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' orbitals added or removed ',(TI(ISTMPO(ISTRIN)-1+I),I=1,ISTMPL(ISTRIN))
  end do
end if
#endif

return

end subroutine CRESTR
