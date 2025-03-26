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

subroutine NSTRSO_MCLR(NEL,NORB1,NORB2,NORB3,NELMN1,NELMX1,NELMN3,NELMX3,IOC,NORB,NSTASO,NOCTYP,NSMST,IOTYP)
! Number of strings per type and symmetry
!
! Jeppe Olsen Winter of 1990

implicit real*8(A-H,O-Z)
dimension IOC(*), NSTASO(NOCTYP,NSMST)

call iCopy(NSMST*NOCTYP,[0],0,NSTASO,1)
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
        IOC(1:IEL1) = [(i,i=1,IEL1)]
        IFRST1 = 0
      else
        call NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1)
        if (NONEW1 == 1) goto 1003
      end if
    end if
    IFRST2 = 1
    IFRST3 = 1
    ! Loop over RAS 2 occupancies
902 continue
    if (IEL2 /= 0) then
      if (IFRST2 == 1) then
        IOC(IEL1+1:IEL1+IEL2) = [(i,i=IORB2F,IORB2F+IEL2-1)]
        IFRST2 = 0
      else
        call NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
        if (NONEW2 == 1) then
          if (IEL1 /= 0) goto 901
          if (IEL1 == 0) goto 1003
        end if
      end if
    end if
    IFRST3 = 1
    ! Loop over RAS 3 occupancies
903 continue
    if (IEL3 /= 0) then
      if (IFRST3 == 1) then
        IOC(IEL1+IEL2+1:IEL1+IEL2+IEL3) = [(i,i=IORB3F,IORB3F+IEL3-1)]
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
    ! Next string has been constructed, Enlist it!
    NSTRIN = NSTRIN+1
    ! Symmetry of string
    ISYM = ISYMST_MCLR(IOC,NEL)
    !      ISYMST_MCLR(STRING,NEL)
    ! occupation type of string
    ITYP = IOCTP2_MCLR(IOC,NEL,IOTYP)
    !      IOCTP2_MCLR(STRING,NEL)

    NSTASO(ITYP,ISYM) = NSTASO(ITYP,ISYM)+1

    if (IEL3 /= 0) goto 903
    if ((IEL3 == 0) .and. (IEL2 /= 0)) goto 902
    if ((IEL3 == 0) .and. (IEL2 == 0) .and. (IEL1 /= 0)) goto 901
1003 continue
  end do
1001 continue
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(NORB)

end subroutine NSTRSO_MCLR
