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
! Copyright (C) 2012, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NEXT_SYM_DISTR_NEW(NSMST,NGRP,KGRP,NGAS,ISYM,ISYM_TOT,IFIRST,NONEW,ISMDFGP,NACTSYM,ISMSCR)
! Giovanni Li Manni. Geneva, February 2012
!
! Obtain next distribution of symmetries with given total
! Symmetry.
!
! Loop over first NGAS-1 spaces are performed, and the symmetry
! of the last space is then fixed by the required total sym

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMST, NGRP, NGAS, KGRP(NGAS), ISYM_TOT, ISMDFGP(NSMST,NGRP), NACTSYM(NGRP)
integer(kind=iwp), intent(out) :: ISYM(NGAS), NONEW
integer(kind=iwp), intent(inout) :: IFIRST, ISMSCR(NGRP)
integer(kind=iwp) :: IGAS, JSYM
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I, ISM
#endif
integer(kind=iwp), external :: ISYMSTR

#ifdef _DEBUGPRINT_
write(u6,*) '***************************'
write(u6,*) 'INPUT IN NEXT_SYM_DISTR_NEW'
write(u6,*) '***************************'
write(u6,*) 'NGRP :',NGRP
write(u6,*) 'KGRP :'
write(u6,'(40I3)') (KGRP(i),i=1,NGAS)
write(u6,*) 'NACTSYM(NGRP) :'
write(u6,'(40I2)') (NACTSYM(i),i=1,NGRP)
write(u6,*) 'ISMDFGP(ISYM,NGRP) :'
do Ism=1,NSMST
  write(u6,'(40I2)') (ISMDFGP(ism,i),i=1,NGRP)
end do
write(u6,*) 'IFIRST',IFIRST
#endif

if (IFIRST == 1) then
  ISMSCR(1:NGAS) = 1
  do IGAS=1,NGAS
    ISYM(IGAS) = ISMDFGP(1,KGRP(IGAS))
  end do
  NONEW = 0
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'Symmetry distribution :'
write(u6,'(40I2)') (ISYM(IGAS),IGAS=1,NGAS)
#endif

do
  if (IFIRST == 0) then
    call NXTDIST(NGRP,NGAS,KGRP,ISMSCR,NACTSYM,NONEW)
    do IGAS=1,NGAS
      ISYM(IGAS) = ISMDFGP(ISMSCR(IGAS),KGRP(IGAS))
    end do
  end if
  IFIRST = 0
  if (NONEW == 0) then
    JSYM = ISYMSTR(ISYM,NGAS)
    if (JSYM == ISYM_TOT) exit
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

end subroutine NEXT_SYM_DISTR_NEW
