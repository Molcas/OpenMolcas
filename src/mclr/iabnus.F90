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

integer function IABNUS(IASTR,NAEL,IAORD,ITPFSA,ISMFSA,NOCTPA,ZA,ISSOA,NSSOA,IBSTR,NBEL,IBORD,ITPFSB,ISMFSB,NOCTPB,ZB,ISSOB,NSSOB, &
                        IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,PSSIGN,IPSFAC,IPRNT)
! A determinant is given by strings IASTR,IBSTR.
! Find number of this determinant
!
! If PSSIGN /= 0, the determinant with higher alpha number is picked
! and phase factor IPSFAC calculated. This corresponds to
! configuration order
!
! Jeppe Olsen

implicit none
integer NAEL
integer IASTR(NAEL), IAORD(*), ITPFSA(*), ISMFSA(*)
integer NOCTPA
integer ZA(*), ISSOA(NOCTPA,*), NSSOA(NOCTPA,*)
integer NBEL
integer IBSTR(NBEL), IBORD(*), ITPFSB(*), ISMFSB(*)
integer NOCTPB
integer ZB(*), ISSOB(NOCTPB,*), NSSOB(NOCTPB,*)
integer IOOS(NOCTPA,NOCTPB,*)
integer NORB, IGENSG
integer ISGNA(*), ISGNB(*)
real*8 PSSIGN
integer IPSFAC, IPRNT
! Local variables
integer NTEST, IANUM, IBNUM, ISGNAB, IASYM, IBSYM, IATP, IBTP, IAREL, IBREL, ISTRNM

NTEST = 000
NTEST = max(NTEST,IPRNT)
if (NTEST > 300) then
  write(6,*) ' >>> IABNUS SPEAKING <<<'
  write(6,*) ' NOCTPA,NOCTPB ',NOCTPA,NOCTPB
  write(6,*) ' ALPHA AND BETA STRING'
  call IWRTMA(IASTR,1,NAEL,1,NAEL)
  call IWRTMA(IBSTR,1,NBEL,1,NBEL)
end if
! Number of alpha- and beta-string
!       ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
IANUM = ISTRNM(IASTR,NORB,NAEL,ZA,IAORD,1)
IBNUM = ISTRNM(IBSTR,NORB,NBEL,ZB,IBORD,1)
if (NTEST >= 10) write(6,*) ' IANUM AND IBNUM ',IANUM,IBNUM

if (IGENSG /= 0) then
  ISGNAB = ISGNA(IANUM)*ISGNB(IBNUM)
else
  ISGNAB = 1
end if
! Symmetries and types
IASYM = ISMFSA(IANUM)
IBSYM = ISMFSB(IBNUM)
!if (NTEST >= 10) write(6,*) ' IASYM IBSYM ',IASYM,IBSYM
IATP = ITPFSA(IANUM)
IBTP = ITPFSB(IBNUM)
!if (NTEST >= 10) write(6,*) ' IATP,IBTP ',IATP,IBTP
IAREL = IANUM-ISSOA(IATP,IASYM)+1
IBREL = IBNUM-ISSOB(IBTP,IBSYM)+1
!if (NTEST >= 10) write(6,*) ' IAREL IBREL ',IAREL,IBREL

if (PSSIGN == 0.0d0) then
  ! Normal determinant ordering
  IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
  IPSFAC = 1
else if (PSSIGN /= 0.0d0) then
  ! Ensure mapping to proper determinant in combination
  if (IANUM >= IBNUM) then
    ! No need for switching around so
    if ((IASYM == IBSYM) .and. (IATP == IBTP)) then
      ! Lower triangular packed, column wise !
      IABNUS = IOOS(IATP,IBTP,IASYM)-1+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-IBREL*(IBREL-1)/2
    else
      IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
    end if
    IPSFAC = 1
  else if (IBNUM > IANUM) then
    ! Switch alpha and beta string around
    if ((IASYM == IBSYM) .and. (IATP == IBTP)) then
      ! Lower triangular packed, column wise !
      IABNUS = IOOS(IBTP,IATP,IBSYM)-1+(IAREL-1)*NSSOB(IBTP,IBSYM)+IBREL-IAREL*(IAREL-1)/2
    else
      IABNUS = IOOS(IBTP,IATP,IBSYM)+(IAREL-1)*NSSOB(IBTP,IBSYM)+IBREL-1
    end if
    IPSFAC = nint(PSSIGN)
  end if

end if

!OLD
!OLD IABNUS = IOOS(IATP,IBTP,IASYM)+(IBREL-1)*NSSOA(IATP,IASYM)+IAREL-1
!if (NTEST > 10) write(6,*) ' IOOS NSSOA ',IOOS(IATP,IBTP,IASYM),NSSOA(IATP,IASYM)

if (NTEST >= 200) then
  write(6,*) ' ALPHA AND BETA STRING'
  call IWRTMA(IASTR,1,NAEL,1,NAEL)
  call IWRTMA(IBSTR,1,NBEL,1,NBEL)
  write(6,*) ' Corresponding determinant number ',IABNUS
end if

end function IABNUS
