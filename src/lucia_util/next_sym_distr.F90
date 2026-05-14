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

!#define _DEBUGPRINT_
subroutine NEXT_SYM_DISTR(NGAS,MNVAL,MXVAL,ISYM,ISYM_TOT,IFIRST,NONEW)
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

use Symmetry_Info, only: Mul
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NGAS, MNVAL(NGAS), MXVAL(NGAS), ISYM_TOT
integer(kind=iwp), intent(inout) :: ISYM(NGAS), IFIRST
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: JSYM
integer(kind=iwp), external :: ISYMSTR

! Symmetry of first NGAS -1 spaces

if (IFIRST == 1) then
  ISYM(1:NGAS-1) = MNVAL(1:NGAS-1)
  NONEW = 0
end if
do
  if (IFIRST == 0) call NXTNUM3(ISYM,NGAS-1,MNVAL,MXVAL,NONEW)
  IFIRST = 0

  ! Symmetry of last space

  if (NONEW == 0) then
    !JSYM = 1
    !do IGAS=1,NGAS-1
    !  JSYM = Mul(JSYM,ISYM(IGAS))
    !end do
    JSYM = ISYMSTR(ISYM,NGAS-1)
    ISYM(NGAS) = Mul(JSYM,ISYM_TOT)

    if ((MNVAL(NGAS) <= ISYM(NGAS)) .and. (MXVAL(NGAS) >= ISYM(NGAS))) exit
  else
    exit
  end if
end do

#ifdef _DEBUGPRINT_
if (NONEW == 1) then
  write(u6,*) ' No new symmetry distributions'
else
  write(u6,*) ' Next symmetry distribution'
  call IWRTMA(ISYM,1,NGAS,1,NGAS)
end if
#endif

end subroutine NEXT_SYM_DISTR
