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

subroutine NSTRSO_MCLR(NEL,NORB1,NORB2,NORB3,NELMN1,NELMX1,NELMN3,NELMX3,IOC,NSTASO,NOCTYP,NSM,IOTYP)
! Number of strings per type and symmetry
!
! Jeppe Olsen Winter of 1990

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NEL, NORB1, NORB2, NORB3, NELMN1, NELMX1, NELMN3, NELMX3, NOCTYP, NSM, IOTYP
integer(kind=iwp), intent(out) :: IOC(NEL), NSTASO(NOCTYP,NSM)
integer(kind=iwp) :: i, IEL1, IEL2, IEL3, IFRST1, IFRST2, IFRST3, IORB1F, IORB1L, IORB2F, IORB2L, IORB3F, IORB3L, ISYM, ITYP, &
                     NONEW1, NONEW2, NONEW3, NSTRIN
integer(kind=iwp), external :: ISYMST_MCLR, IOCTP2_MCLR

NSTASO(:,:) = 0
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
          ! Next string has been constructed, Enlist it!
          NSTRIN = NSTRIN+1
          ! Symmetry of string
          ISYM = ISYMST_MCLR(IOC,NEL)
          !      ISYMST_MCLR(STRING,NEL)
          ! occupation type of string
          ITYP = IOCTP2_MCLR(IOC,NEL,IOTYP)
          !      IOCTP2_MCLR(STRING,NEL)

          NSTASO(ITYP,ISYM) = NSTASO(ITYP,ISYM)+1

          if (IEL3 == 0) exit RAS3occ
        end do RAS3occ
        if ((IEL3 /= 0) .or. (IEL2 == 0)) exit RAS2occ
      end do RAS2occ
      if ((IEL3 /= 0) .or. (IEL2 /= 0) .or. (IEL1 == 0)) exit RAS1occ
    end do RAS1occ
  end do inner
end do outer

end subroutine NSTRSO_MCLR
