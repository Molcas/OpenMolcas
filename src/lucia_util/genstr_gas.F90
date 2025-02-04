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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************

subroutine GENSTR_GAS(NEL,NELMN1,NELMX1,NELMN3,NELMX3,ISTASO,IGRP,NOCTYP,NSMST,Z,LSTASO,IREORD,STRING,IOC,IPRNT)
! Generate strings consisting of  NEL electrons fulfilling
!   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
!   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
!
! In the present version the strings are directly ordered into
! symmetry and occupation type .
!
! Jeppe Olsen Winter of 1990
!
! Special GAS version, Winter of 94 All strings of group IGRP
!
! ========
! Output :
! ========
! STRING(IEL,ISTRIN) : Occupation of strings.
! IREORD             : Reordering array going from lexical
!                      order to symmetry and occupation type order.

use lucia_data, only: NACOB, NORB1, NORB2, NORB3
use Definitions, only: u6

implicit none
integer NEL, NELMN1, NELMX1, NELMN3, NELMX3, IGRP, NOCTYP, NSMST, IPRNT
! Input
integer ISTASO(NSMST,*)
! Orbinp
integer Z(NACOB,NEL)
! Output
integer STRING(NEL,*), IREORD(*)
! Scratch arrays
integer IOC(*), LSTASO(NOCTYP,NSMST)
integer NTEST0, NTEST, NSTRIN, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, IEL1, IEL3, IEL2, IFRST1, NONEW1, IFRST2, IFRST3, &
        NONEW2, NONEW3, ISYM, ITYP, LEXCI, LACTU, NPR, ISTRIN, LSTRIN, KSTRIN, IEL, ISTRNM, ISYMST

call ISETVC(LSTASO,0,NOCTYP*NSMST)
NTEST0 = 0
NTEST = max(NTEST0,IPRNT)
if (NTEST >= 10) then
  write(u6,*) ' ==============='
  write(u6,*) ' GENSTR speaking'
  write(u6,*) ' ==============='
end if

NSTRIN = 0
IORB1F = 1
IORB1L = IORB1F+NORB1-1
IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1
IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1
! Loop over possible partitionings between RAS1,RAS2,RAS3
do IEL1=NELMX1,NELMN1,-1
  do IEL3=NELMN3,NELMX3,1
    if (IEL1 > NORB1) goto 1001
    if (IEL3 > NORB3) goto 1003
    IEL2 = NEL-IEL1-IEL3
    if ((IEL2 < 0) .or. (IEL2 > NORB2)) goto 1003
    IFRST1 = 1
    ! Loop over RAS 1 occupancies
901 continue
    if (IEL1 /= 0) then
      if (IFRST1 == 1) then
        call ISTVC2(IOC(1),0,1,IEL1)
        IFRST1 = 0
      else
        call NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
        if (NONEW1 == 1) goto 1003
      end if
    end if
    if (NTEST >= 500) then
      write(u6,*) ' RAS 1 string'
      call IWRTMA(IOC,1,IEL1,1,IEL1)
    end if
    IFRST2 = 1
    IFRST3 = 1
    ! Loop over RAS 2 occupancies
902 continue
    if (IEL2 /= 0) then
      if (IFRST2 == 1) then
        call ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
        IFRST2 = 0
      else
        call NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
        if (NONEW2 == 1) then
          if (IEL1 /= 0) goto 901
          if (IEL1 == 0) goto 1003
        end if
      end if
    end if
    if (NTEST >= 500) then
      write(u6,*) ' RAS 1 2 string'
      call IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
    end if
    IFRST3 = 1
    ! Loop over RAS 3 occupancies
903 continue
    if (IEL3 /= 0) then
      if (IFRST3 == 1) then
        call ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
        IFRST3 = 0
      else
        call NXTORD(IOC(IEL1+IEL2+1),IEL3,IORB3F,IORB3L,NONEW3)
        if (NONEW3 == 1) then
          if (IEL2 /= 0) goto 902
          if (IEL1 /= 0) goto 901
          goto 1003
        end if
      end if
    end if
    if (NTEST >= 500) then
      write(u6,*) ' RAS 1 2 3 string'
      call IWRTMA(IOC,1,NEL,1,NEL)
    end if
    ! Next string has been constructed, enlist it!
    NSTRIN = NSTRIN+1
    ! Symmetry
    !ISYM = ISYMST(STRING,NEL)
    ISYM = ISYMST(IOC,NEL)
    ! Occupation type
    !ITYP = IOCTP2(IOC,NEL,IOTYP)
    ITYP = 1

    if (ITYP /= 0) then
      LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+1
      !ISTRNM(IOCC,NACTOB,NEL,Z,NEWORD,IREORD)
      LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
      LACTU = ISTASO(ISYM,IGRP)-1+LSTASO(ITYP,ISYM)
      IREORD(LEXCI) = LACTU
      if (NTEST > 10) write(u6,*) ' LEXCI,LACTU',LEXCI,LACTU
      if (NEL > 0) call ICOPVE(IOC,STRING(1,LACTU),NEL)
    end if

    if (IEL3 /= 0) goto 903
    if ((IEL3 == 0) .and. (IEL2 /= 0)) goto 902
    if ((IEL3 == 0) .and. (IEL2 == 0) .and. (IEL1 /= 0)) goto 901
1003 continue
  end do
1001 continue
end do

if (NTEST >= 1) write(u6,*) ' Number of strings generated   ',NSTRIN
if (NTEST >= 10) then
  if (NTEST >= 100) then
    NPR = NSTRIN
  else
    NPR = min(NSTRIN,50)
  end if
  write(u6,*) ' Strings generated'
  write(u6,*) ' ================='
  ISTRIN = 0
  do ISYM=1,NSMST
    do ITYP=1,NOCTYP
      LSTRIN = min(LSTASO(ITYP,ISYM),NPR-ISTRIN)
      if (LSTRIN > 0) then
        write(u6,*) ' Strings of type and symmetry ',ITYP,ISYM
        do KSTRIN=1,LSTRIN
          ISTRIN = ISTRIN+1
          write(u6,'(2X,I4,8X,(10I5))') ISTRIN,(STRING(IEL,ISTRIN),IEL=1,NEL)
        end do
      end if
    end do
  end do

  write(u6,*) ' Array giving actual place from lexical place'
  write(u6,*) ' ============================================'
  call IWRTMA(IREORD,1,NPR,1,NPR)
end if

end subroutine GENSTR_GAS
