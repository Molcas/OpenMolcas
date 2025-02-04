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

subroutine NSTRSO_GAS(NEL,NORB1,NORB2,NORB3,NELMN1,NELMX1,NELMN3,NELMX3,IOC,NORB,NSTASO,ISTASO,NOCTYP,NSMST,IOTYP,IPRNT)
! Number of strings per symmetry for group IOTYP
!
! Gas version, no check of type : set to 1
!
! Jeppe Olsen Winter of 1994

use Definitions, only: u6

implicit real*8(A-H,O-Z)
!dimension IOC(*), NSTASO(NOCTYP,NSMST)
dimension IOC(*), NSTASO(NSMST,*), ISTASO(NSMST,*)

call ISETVC(NSTASO(1,IOTYP),0,NSMST)
NTEST0 = 0
NTEST = max(IPRNT,NTEST0)
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
    ! Next string has been constructed, enlist it !
    NSTRIN = NSTRIN+1
    ! Symmetry of string
    ISYM = ISYMST(IOC,NEL)
    !      ISYMST(STRING,NEL)
    ! occupation type of string
    !OLD ITYP = IOCTP2(IOC,NEL,IOTYP)
    !           IOCTP2(STRING,NEL)

    NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+1

    if (IEL3 /= 0) goto 903
    if ((IEL3 == 0) .and. (IEL2 /= 0)) goto 902
    if ((IEL3 == 0) .and. (IEL2 == 0) .and. (IEL1 /= 0)) goto 901
1003 continue
  end do
1001 continue
end do

! The corresponding offset

do ISM=1,NSMST
  if (ISM == 1) then
    ISTASO(ISM,IOTYP) = 1
  else
    ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
  end if
end do

if (NTEST >= 5) write(u6,*) ' Number of strings generated   ',NSTRIN
if (NTEST >= 10) then
  write(u6,*)
  write(u6,*) ' Number of strings per sym for group = ',IOTYP
  write(u6,*) '================================================'
  call IWRTMA(NSTASO(1,IOTYP),1,NSMST,1,NSMST)
  write(u6,*) ' Offset for given symmetry for group = ',IOTYP
  write(u6,*) '================================================'
  call IWRTMA(ISTASO(1,IOTYP),1,NSMST,1,NSMST)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(NORB)
  call Unused_integer(NOCTYP)
end if

end subroutine NSTRSO_GAS
