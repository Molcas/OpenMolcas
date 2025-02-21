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
subroutine CRESTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,Z,NEWORD,TI,TTO,NACOB)
! A type of strings containing NEL electrons are given
! set up all possible ways of adding an electron to this type of strings
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
!          NEL + 1 electrons
! NEWORD : Reordering array for N+1 strings
! LSGSTR : /= 0 => Include sign arrays ISGSTI,ISGSTO of strings
! ISGSTI : Sign array for NEL   strings
! ISGSTO : Sign array for NEL+1 strings
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

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSTINI, NEL, STRING(NEL,NSTINI), NSTINO, NORB, IORBOF, Z(NORB,NEL+1), NEWORD(NSTINO), NACOB
integer(kind=iwp), intent(out) :: TI(NORB,NSTINI), TTO(NORB,NSTINI)
integer(kind=iwp) :: IEL, IIISGN, IORB, IPLACE, ISTRIN, JSTRIN, STRIN2(500)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I, NPR
#endif
integer(kind=iwp), external :: ISTRNM

#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' CRESTR speaking'
write(u6,*) ' ==============='
write(u6,*)
write(u6,*) ' Number of input electrons ',NEL
!write(u6,*) ' Reorder array NEWORD'
!call IWRTMA(NEWORD,1,NSTINO,1,NSTINO)
#endif

do ISTRIN=1,NSTINI
  do IORB=IORBOF,IORBOF-1+NORB

    IPLACE = 0

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
    !write(u6,*) ' updated string (STRIN2)'
    !call iwrtma(STRIN2,1,NEL+1,1,NEL+1)
    JSTRIN = ISTRNM(STRIN2,NACOB,NEL+1,Z,NEWORD,1)
    !write(u6,*) ' corresponding number ',JSTRIN

    TTO(IORB-IORBOF+1,ISTRIN) = JSTRIN
    IIISGN = (-1)**(IPLACE-1)
    !if (LSGSTR /= 0) IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
    if (IIISGN == -1) TTO(IORB-IORBOF+1,ISTRIN) = -TTO(IORB-IORBOF+1,ISTRIN)
    TI(IORB-IORBOF+1,ISTRIN) = IORB

  end do

end do

#ifdef _DEBUGPRINT_
NPR = min(NSTINI,60)
write(u6,*) ' Output from CRESTR :'
write(u6,*) '==================='

write(u6,*)
write(u6,*) ' Strings with an electron added'
do ISTRIN=1,NPR
  write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' New strings.. ',(TTO(I,ISTRIN),I=1,NORB)
end do
do ISTRIN=1,NPR
  write(u6,'(2X,A,I4,A,/,(10I5))') 'String..',ISTRIN,' orbitals added or removed ',(TI(I,ISTRIN),I=1,NORB)
end do
#endif

end subroutine CRESTR_GAS
