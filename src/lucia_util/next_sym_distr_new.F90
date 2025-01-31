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

subroutine NEXT_SYM_DISTR_NEW(NSMST,NGRP,KGRP,NGAS,ISYM,ISYM_TOT,IFIRST,NONEW,ISMDFGP,NACTSYM,ISMSCR)
! Giovanni Li Manni. Geneva, February 2012
!
! Obtain next distribution of symmetries with given total
! Symmetry.
!
! Loop over first NGAS-1 spaces are performed, and the symmetry
! of the last space is then fixed by the required total sym

implicit real*8(A-H,O-Z)
! Input
integer ISMDFGP(NSMST,NGRP), NACTSYM(NGRP), ISMSCR(NGRP)
integer KGRP(NGAS)
! Input and output
integer ISYM(NGAS)

NTEST = 0

if (NTEST >= 100) then
  write(6,*) '***************************'
  write(6,*) 'INPUT IN NEXT_SYM_DISTR_NEW'
  write(6,*) '***************************'
  write(6,*) 'NGRP :',NGRP
  write(6,*) 'KGRP :'
  write(6,'(40I3)') (KGRP(i),i=1,NGAS)
  write(6,*) 'NACTSYM(NGRP) :'
  write(6,'(40I2)') (NACTSYM(i),i=1,NGRP)
  write(6,*) 'ISMDFGP(ISYM,NGRP) :'
  do Ism=1,NSMST
    write(6,'(40I2)') (ISMDFGP(ism,i),i=1,NGRP)
  end do
  write(6,*) 'IFIRST',IFIRST
end if

if (IFIRST == 1) then
  do IGAS=1,NGAS
    ISMSCR(IGAS) = 1
    ISYM(IGAS) = ISMDFGP(1,KGRP(IGAS))
  end do
  NONEW = 0
end if

if (NTEST >= 100) then
  write(6,*) 'Symmetry distribution :'
  write(6,'(40I2)') (ISYM(IGAS),IGAS=1,NGAS)
end if

1001 continue
if (IFIRST == 0) then
  call NXTDIST(NSMST,NGRP,NGAS,KGRP,ISMDFGP,ISMSCR,NACTSYM,NONEW)
  do IGAS=1,NGAS
    ISYM(IGAS) = ISMDFGP(ISMSCR(IGAS),KGRP(IGAS))
  end do
end if
IFIRST = 0
if (NONEW == 0) then
  JSYM = ISYMSTR(ISYM,NGAS)
  if (JSYM /= ISYM_TOT) goto 1001
end if

if (NTEST >= 100) then
  if (NONEW == 1) then
    write(6,*) ' No new symmetry distributions'
  else
    write(6,*) ' Next symmetry distribution'
    call IWRTMA(ISYM,1,NGAS,1,NGAS)
  end if
end if

end subroutine NEXT_SYM_DISTR_NEW
