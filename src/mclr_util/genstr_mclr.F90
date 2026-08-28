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

subroutine GENSTR_MCLR(NEL,NELMN1,NELMX1,NELMN3,NELMX3,ISTASO,NOCTYP,NSM,Z,LSTASO,IREORD,STRING,IOC,IOTYP)
! Generate strings consisting of  NEL electrons fulfilling
!   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
!   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
!
! In the present version the strings are directly ordered into
! symmetry and occupation type.
!
! Jeppe Olsen Winter of 1990
! ========
! Output :
! ========
! STRING(IEL,ISTRIN) : Occupation of strings.
! IREORD             : Reordering array going from lexical
!                      order to symmetry and occupation type order.

use MCLR_Data, only: NACOB, NORB1, NORB2, NORB3
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NEL, NELMN1, NELMX1, NELMN3, NELMX3, NOCTYP, NSM, ISTASO(NOCTYP,NSM), Z(NACOB,NEL), IOTYP
integer(kind=iwp), intent(out) :: LSTASO(NOCTYP,NSM)
integer(kind=iwp), intent(inout) :: IREORD(*)
integer(kind=iwp), intent(_OUT_) :: STRING(NEL,*), IOC(*)
integer(kind=iwp) :: i, IEL1, IEL2, IEL3, IFRST1, IFRST2, IFRST3, IOCTP2_MCLR, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, &
                     ISTRNM, ISYM, ISYMST_MCLR, ITYP, LACTU, LEXCI, NONEW1, NONEW2, NONEW3, NSTRIN
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IEL, ISTRIN, KSTRIN, LSTRIN
#endif

LSTASO(:,:) = 0
NSTRIN = 0
IORB1F = 1
IORB1L = IORB1F+NORB1-1
IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1
IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1
! Loop over possible partitionings between RAS1,RAS2,RAS3
outer: do IEL1=NELMX1,NELMN1,-1
  inner: do IEL3=NELMN3,NELMX3,1
    if (IEL1 > NORB1) cycle outer
    if (IEL3 > NORB3) cycle inner
    IEL2 = NEL-IEL1-IEL3
    if ((IEL2 < 0) .or. (IEL2 > NORB2)) cycle inner
    IFRST1 = 1
    ! Loop over RAS 1 occupancies
    RAS1occ: do
      if (IEL1 /= 0) then
        if (IFRST1 == 1) then
          IOC(1:IEL1) = [(i,i=1,IEL1)]
          IFRST1 = 0
        else
          call NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
          if (NONEW1 == 1) cycle inner
        end if
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) ' RAS 1 string'
      call IWRTMA(IOC,1,IEL1,1,IEL1)
#     endif
      IFRST2 = 1
      IFRST3 = 1
      ! Loop over RAS 2 occupancies
      RAS2occ: do
        if (IEL2 /= 0) then
          if (IFRST2 == 1) then
            IOC(IEL1+1:IEL1+IEL2) = [(i,i=IORB2F,IORB2F+IEL2-1)]
            IFRST2 = 0
          else
            call NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
            if (NONEW2 == 1) then
              if (IEL1 /= 0) cycle RAS1occ
              if (IEL1 == 0) cycle inner
            end if
          end if
        end if
#       ifdef _DEBUGPRINT_
        write(u6,*) ' RAS 1 2 string'
        call IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
#       endif
        IFRST3 = 1
        ! Loop over RAS 3 occupancies
        RAS3occ: do
          if (IEL3 /= 0) then
            if (IFRST3 == 1) then
              IOC(IEL1+IEL2+1:IEL1+IEL2+IEL3) = [(i,i=IORB3F,IORB3F+IEL3-1)]
              IFRST3 = 0
            else
              call NXTORD(IOC(IEL1+IEL2+1),IEL3,IORB3F,IORB3L,NONEW3)
              if (NONEW3 == 1) then
                if (IEL2 /= 0) cycle RAS2occ
                if (IEL1 /= 0) cycle RAS1occ
                cycle inner
              end if
            end if
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) ' RAS 1 2 3 string'
          call IWRTMA(IOC,1,NEL,1,NEL)
#         endif
          ! Next string has been constructed, Enlist it!
          NSTRIN = NSTRIN+1
          ! Symmetry
          ISYM = ISYMST_MCLR(IOC,NEL)
          ! Occupation type
          ITYP = IOCTP2_MCLR(IOC,NEL,IOTYP)

          if (ITYP /= 0) then
            LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+1
            LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
            LACTU = ISTASO(ITYP,ISYM)-1+LSTASO(ITYP,ISYM)
            IREORD(LEXCI) = LACTU
            STRING(:,LACTU) = IOC(1:NEL)
          end if

          if (IEL3 == 0) exit RAS3occ
        end do RAS3occ
        if ((IEL3 /= 0) .or. (IEL2 == 0)) exit RAS2occ
      end do RAS2occ
      if ((IEL3 /= 0) .or. (IEL2 /= 0) .or. (IEL1 == 0)) exit RAS1occ
    end do RAS1occ
  end do inner
end do outer

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings generated ',NSTRIN
write(u6,*) ' Strings generated'
write(u6,*) ' =================='
ISTRIN = 0
do ISYM=1,NSM
  do ITYP=1,NOCTYP
    LSTRIN = min(LSTASO(ITYP,ISYM),NSTRIN-ISTRIN)
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
call IWRTMA(IREORD,1,NSTRIN,1,NSTRIN)
#endif

end subroutine GENSTR_MCLR
