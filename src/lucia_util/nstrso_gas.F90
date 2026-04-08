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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NSTRSO_GAS(NEL,NORB1,NORB2,NORB3,NELMN1,NELMX1,NELMN3,NELMX3,IOC,NSTASO,ISTASO,NSMST,IOTYP)
! Number of strings per symmetry for group IOTYP
!
! Gas version, no check of type : set to 1
!
! Jeppe Olsen Winter of 1994

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NEL, NORB1, NORB2, NORB3, NELMN1, NELMX1, NELMN3, NELMX3, NSMST, IOTYP
integer(kind=iwp), intent(out) :: IOC(NEL)
integer(kind=iwp), intent(inout) :: NSTASO(NSMST,IOTYP), ISTASO(NSMST,IOTYP)
integer(kind=iwp) :: i, IEL1, IEL2, IEL3, IFRST1, IFRST2, IFRST3, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, ISM, ISYM, &
                     NONEW1, NONEW2, NONEW3, NSTRIN
integer(kind=iwp), external :: ISYMST

NSTASO(:,IOTYP) = 0
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
          call NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
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
            call NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
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
              call NXTORD(IOC(IEL1+IEL2+1),IEL3,IORB3F,IORB3L,NONEW3)
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
          ! Next string has been constructed, enlist it !
          NSTRIN = NSTRIN+1
          ! Symmetry of string
          ISYM = ISYMST(IOC,NEL)
          !      ISYMST(STRING,NEL)
          ! occupation type of string
          !OLD ITYP = IOCTP2(IOC,NEL,IOTYP)
          !           IOCTP2(STRING,NEL)

          NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+1

          if (IEL3 == 0) exit ras3
        end do ras3
        if ((IEL3 /= 0) .or. (IEL2 == 0)) exit ras2
      end do ras2
      if ((IEL3 /= 0) .or. (IEL2 /= 0) .or. (IEL1 == 0)) exit ras1
    end do ras1
  end do outer
end do

! The corresponding offset

do ISM=1,NSMST
  if (ISM == 1) then
    ISTASO(ISM,IOTYP) = 1
  else
    ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings generated   ',NSTRIN
write(u6,*)
write(u6,*) ' Number of strings per sym for group = ',IOTYP
write(u6,*) '================================================'
call IWRTMA(NSTASO(:,IOTYP),1,NSMST,1,NSMST)
write(u6,*) ' Offset for given symmetry for group = ',IOTYP
write(u6,*) '================================================'
call IWRTMA(ISTASO(:,IOTYP),1,NSMST,1,NSMST)
#endif

end subroutine NSTRSO_GAS
