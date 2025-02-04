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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************

subroutine NEXT_SYM_DISTR(NGAS,MINVAL,MAXVAL,ISYM,ISYM_TOT,IFIRST,NONEW)
! Obtain next distribution of symmetries with given total
! Symmetry.
!
! Loop over first NGAS-1 spaces are performed, and the symmetry
! of the last space is then fixed by the required total sym
!
! Jeppe Olsen, Sept 97
! Obtain next distribution of symmetries with given total
! Symmetry.
!
! Loop over first NGAS-1 spaces are performed, and the symmetry
! of the last space is then fixed by the required total sym
!
! Jeppe Olsen, Sept 97

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
dimension minval(NGAS), maxval(NGAS)
! Input and output
dimension ISYM(NGAS)

! Symmetry of first NGAS -1 spaces

if (IFIRST == 1) then
  do IGAS=1,NGAS-1
    ISYM(IGAS) = minval(IGAS)
  end do
  NONEW = 0
end if
1001 continue
if (IFIRST == 0) call NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
IFIRST = 0

! Symmetry of last space

if (NONEW == 0) then
  !JSYM = 1
  !do IGAS=1,NGAS-1
  !  call SYMCOM(3,0,JSYM,ISYM(IGAS),KSYM)
  !  JSYM = KSYM
  !end do
  JSYM = ISYMSTR(ISYM,NGAS-1)
  call SYMCOM(2,0,JSYM,ISYM(NGAS),ISYM_TOT)

  if ((minval(NGAS) > ISYM(NGAS)) .or. (maxval(NGAS) < ISYM(NGAS))) goto 1001
end if

NTEST = 0
if (NTEST >= 100) then
  if (NONEW == 1) then
    write(u6,*) ' No new symmetry distributions'
  else
    write(u6,*) ' Next symmetry distribution'
    call IWRTMA(ISYM,1,NGAS,1,NGAS)
  end if
end if

end subroutine NEXT_SYM_DISTR
