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

!#define _DEBUGPRINT_
subroutine GENSTR_GAS(NEL,NELMN1,NELMX1,NELMN3,NELMX3,ISTASO,IGRP,NOCTYP,NSMST,Z,LSTASO,IREORD,STRING,IOC)
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
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NEL, NELMN1, NELMX1, NELMN3, NELMX3, NSMST, ISTASO(NSMST,*), IGRP, NOCTYP, Z(NACOB,NEL)
integer(kind=iwp), intent(out) :: LSTASO(NOCTYP,NSMST), IOC(NEL)
integer(kind=iwp), intent(inout) :: IREORD(*)
integer(kind=iwp), intent(_OUT_) :: STRING(NEL,*)
integer(kind=iwp) :: i, IEL1, IEL2, IEL3, IFRST1, IFRST2, IFRST3, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, ISTRNM, ISYM, &
                     ISYMST, ITYP, LACTU, LEXCI, NONEW1, NONEW2, NONEW3, NSTRIN
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IEL, ISTRIN, KSTRIN, LSTRIN, NPR
#endif

LSTASO(:,:) = 0
#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' GENSTR speaking'
write(u6,*) ' ==============='
#endif

NSTRIN = 0
IORB1F = 1
IORB1L = IORB1F+NORB1-1
IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1
IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1
! Loop over possible partitionings between RAS1,RAS2,RAS3
do IEL1=NELMX1,NELMN1,-1
  outer: do IEL3=NELMN3,NELMX3,1
    if (IEL1 > NORB1) exit outer
    if (IEL3 > NORB3) cycle outer
    IEL2 = NEL-IEL1-IEL3
    if ((IEL2 < 0) .or. (IEL2 > NORB2)) cycle outer
    IFRST1 = 1
    ! Loop over RAS 1 occupancies
    ras1: do
      if (IEL1 /= 0) then
        if (IFRST1 == 1) then
          IOC(1:IEL1) = [(i,i=1,IEL1)]
          IFRST1 = 0
        else
          call NXTORD(IOC(1:IEL1),IEL1,IORB1F,IORB1L,NONEW1)
          if (NONEW1 == 1) cycle outer
        end if
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) ' RAS 1 string'
      call IWRTMA(IOC,1,IEL1,1,IEL1)
#     endif
      IFRST2 = 1
      IFRST3 = 1
      ! Loop over RAS 2 occupancies
      ras2: do
        if (IEL2 /= 0) then
          if (IFRST2 == 1) then
            IOC(IEL1+1:IEL1+IEL2) = [(i,i=IORB2F,IORB2F+IEL2-1)]
            IFRST2 = 0
          else
            call NXTORD(IOC(IEL1+1:IEL1+IEL2),IEL2,IORB2F,IORB2L,NONEW2)
            if (NONEW2 == 1) then
              if (IEL1 /= 0) cycle ras1
              if (IEL1 == 0) cycle outer
            end if
          end if
        end if
#       ifdef _DEBUGPRINT_
        write(u6,*) ' RAS 1 2 string'
        call IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
#       endif
        IFRST3 = 1
        ! Loop over RAS 3 occupancies
        ras3: do
          if (IEL3 /= 0) then
            if (IFRST3 == 1) then
              IOC(IEL1+IEL2+1:IEL1+IEL2+IEL3) = [(i,i=IORB3F,IORB3F+IEL3-1)]
              IFRST3 = 0
            else
              call NXTORD(IOC(IEL1+IEL2+1:IEL1+IEL2+IEL3),IEL3,IORB3F,IORB3L,NONEW3)
              if (NONEW3 == 1) then
                if (IEL2 /= 0) cycle ras2
                if (IEL1 /= 0) cycle ras1
                cycle outer
              end if
            end if
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) ' RAS 1 2 3 string'
          call IWRTMA(IOC,1,NEL,1,NEL)
#         endif
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
            !       ISTRNM(IOCC,NACTOB,NEL,Z,NEWORD,IREORD)
            LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
            LACTU = ISTASO(ISYM,IGRP)-1+LSTASO(ITYP,ISYM)
            IREORD(LEXCI) = LACTU
#           ifdef _DEBUGPRINT_
            write(u6,*) ' LEXCI,LACTU',LEXCI,LACTU
#           endif
            STRING(:,LACTU) = IOC(1:NEL)
          end if

          if (IEL3 == 0) exit ras3
        end do ras3
        if ((IEL3 /= 0) .or. (IEL2 == 0)) exit ras2
      end do ras2
      if ((IEL3 /= 0) .or. (IEL2 /= 0) .or. (IEL1 == 0)) exit ras1
    end do ras1
  end do outer
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings generated   ',NSTRIN
NPR = NSTRIN
!NPR = min(NSTRIN,50)
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
#endif

end subroutine GENSTR_GAS
