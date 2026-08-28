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

function IABNUS(IASTR,NAEL,IAORD,ITPFSA,ISMFSA,NOCTPA,ZA,ISSOA,NSSOA,IBSTR,NBEL,IBORD,ITPFSB,ISMFSB,NOCTPB,ZB,ISSOB,NSSOB,IOOS, &
                NORB,IGENSG,ISGNA,ISGNB,ISGNAB,PSSIGN,IPSFAC)
! A determinant is given by strings IASTR,IBSTR.
! Find number of this determinant
!
! If PSSIGN /= 0, the determinant with higher alpha number is picked
! and phase factor IPSFAC calculated. This corresponds to
! configuration order
!
! Jeppe Olsen

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IABNUS
integer(kind=iwp), intent(in) :: NAEL, IASTR(NAEL), IAORD(*), ITPFSA(*), ISMFSA(*), NOCTPA, ZA(*), ISSOA(NOCTPA,*), &
                                 NSSOA(NOCTPA,*), NBEL, IBSTR(NBEL), IBORD(*), ITPFSB(*), ISMFSB(*), NOCTPB, ZB(*), &
                                 ISSOB(NOCTPB,*), NSSOB(NOCTPB,*), IOOS(NOCTPA,NOCTPB,*), NORB, IGENSG, ISGNA(*), ISGNB(*)
integer(kind=iwp), intent(out) :: ISGNAB, IPSFAC
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IANUM, IAREL, IASYM, IATP, IBNUM, IBREL, IBSYM, IBTP, ISTRNM

#ifdef _DEBUGPRINT_
write(u6,*) ' >>> IABNUS SPEAKING <<<'
write(u6,*) ' NOCTPA,NOCTPB ',NOCTPA,NOCTPB
write(u6,*) ' ALPHA AND BETA STRING'
call IWRTMA(IASTR,1,NAEL,1,NAEL)
call IWRTMA(IBSTR,1,NBEL,1,NBEL)
#endif
! Number of alpha- and beta-string
!       ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
IANUM = ISTRNM(IASTR,NORB,NAEL,ZA,IAORD,1)
IBNUM = ISTRNM(IBSTR,NORB,NBEL,ZB,IBORD,1)
#ifdef _DEBUGPRINT_
write(u6,*) ' IANUM AND IBNUM ',IANUM,IBNUM
#endif

if (IGENSG /= 0) then
  ISGNAB = ISGNA(IANUM)*ISGNB(IBNUM)
else
  ISGNAB = 1
end if
! Symmetries and types
IASYM = ISMFSA(IANUM)
IBSYM = ISMFSB(IBNUM)
!#ifdef _DEBUGPRINT_
!write(u6,*) ' IASYM IBSYM ',IASYM,IBSYM
!#endif
IATP = ITPFSA(IANUM)
IBTP = ITPFSB(IBNUM)
!#ifdef _DEBUGPRINT_
!write(u6,*) ' IATP,IBTP ',IATP,IBTP
!#endif
IAREL = IANUM-ISSOA(IATP,IASYM)+1
IBREL = IBNUM-ISSOB(IBTP,IBSYM)+1
!#ifdef _DEBUGPRINT_
!write(u6,*) ' IAREL IBREL ',IAREL,IBREL
!#endif

if (PSSIGN == Zero) then
  ! Normal determinant ordering
  IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
  IPSFAC = 1
else if (PSSIGN /= Zero) then
  ! Ensure mapping to proper determinant in combination
  if (IANUM >= IBNUM) then
    ! No need for switching around so
    if ((IASYM == IBSYM) .and. (IATP == IBTP)) then
      ! Lower triangular packed, column wise !
      IABNUS = IOOS(IATP,IBTP,IASYM)-1+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-nTri_Elem(IBREL-1)
    else
      IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
    end if
    IPSFAC = 1
  else if (IBNUM > IANUM) then
    ! Switch alpha and beta string around
    if ((IASYM == IBSYM) .and. (IATP == IBTP)) then
      ! Lower triangular packed, column wise !
      IABNUS = IOOS(IBTP,IATP,IBSYM)-1+(IAREL-1)*NSSOB(IBTP,IBSYM)+IBREL-nTri_Elem(IAREL-1)
    else
      IABNUS = IOOS(IBTP,IATP,IBSYM)+(IAREL-1)*NSSOB(IBTP,IBSYM)+IBREL-1
    end if
    IPSFAC = nint(PSSIGN)
  end if

end if

!OLD
!OLD IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
!#ifdef _DEBUGPRINT_
!write(u6,*) ' IOOS NSSOA ',IOOS(IATP,IBTP,IASYM),NSSOA(IATP,IASYM)
!#endif

#ifdef _DEBUGPRINT_
write(u6,*) ' ALPHA AND BETA STRING'
call IWRTMA(IASTR,1,NAEL,1,NAEL)
call IWRTMA(IBSTR,1,NBEL,1,NBEL)
write(u6,*) ' Corresponding determinant number ',IABNUS
#endif

end function IABNUS
