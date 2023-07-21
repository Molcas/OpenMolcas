!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CHO_PRTRED(IOPT)
!
! Purpose: print information about reduced set.

use ChoSwp, only: nnBstRSh

implicit real*8(a-h,o-z)
#include "choorb.fh"
#include "cholesky.fh"
real*8 XBAS(8), XXBAS(8)
integer NSHP(2)
logical CONTRIB(2)
! Statement function
MULD2H(I,J) = ieor(I-1,J-1)+1

do ISYM=1,NSYM
  XBAS(ISYM) = dble(NBAS(ISYM))
end do

XXBAST = 0.0d0
do ISYM=1,NSYM
  XXBAS(ISYM) = 0.0d0
  do ISYMB=1,NSYM
    ISYMA = MULD2H(ISYMB,ISYM)
    if (ISYMA == ISYMB) then
      XXBAS(ISYM) = XXBAS(ISYM)+XBAS(ISYMA)*(XBAS(ISYMA)+1.0d0)/2.0d0
    else if (ISYMA > ISYMB) then
      XXBAS(ISYM) = XXBAS(ISYM)+XBAS(ISYMA)*XBAS(ISYMB)
    end if
  end do
  XXBAST = XXBAST+XXBAS(ISYM)
end do

if (IOPT == 1) then
  NRED = 1
else
  NRED = 2
end if

do IRED=1,NRED
  NSHP(IRED) = 0
  do ISHLAB=1,NNSHL
    CONTRIB(IRED) = .false.
    ISYM = 0
    do while ((.not. CONTRIB(IRED)) .and. (ISYM < NSYM))
      ISYM = ISYM+1
      if (NNBSTRSH(ISYM,ISHLAB,IRED) > 0) CONTRIB(IRED) = .true.
    end do
    if (CONTRIB(IRED)) NSHP(IRED) = NSHP(IRED)+1
  end do
end do

call CHO_HEAD('Reduced Set Information','=',80,LUPRI)

if (NNSHL_TOT == 0) then
  PCT1 = 9.9d9
else
  PCT1 = 1.0d2*dble(NSHP(1))/dble(NNSHL_TOT)
end if
if (IOPT == 1) then
  write(LUPRI,'(/,A,/,A)') 'Sym.          Full   First Red. Set','-----------------------------------'
  do ISYM=1,NSYM
    write(LUPRI,'(I3,3X,F12.1,7X,I10)') ISYM,XXBAS(ISYM),NNBSTR(ISYM,1)
  end do
  write(LUPRI,'(A)') '-----------------------------------'
  write(LUPRI,'(A,F12.1,7X,I10)') 'Total:',XXBAST,NNBSTRT(1)
  write(LUPRI,'(A)') '-----------------------------------'
  write(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)') 'First Reduced Set:',NSHP(1),' of',NNSHL_TOT,' shell pairs contribute (',PCT1,'%)'
else
  if (NNSHL_TOT == 0) then
    PCT2 = 9.9d9
  else
    PCT2 = 1.0d2*dble(NSHP(2))/dble(NNSHL_TOT)
  end if
  write(LUPRI,'(/,A,/,A,/,A)') '                          Reduced Set','Sym.          Full      First    Current', &
                               '----------------------------------------'
  do ISYM=1,NSYM
    write(LUPRI,'(I3,3X,F12.1,1X,I10,1X,I10)') ISYM,XXBAS(ISYM),NNBSTR(ISYM,1),NNBSTR(ISYM,2)
  end do
  write(LUPRI,'(A)') '----------------------------------------'
  write(LUPRI,'(A,F12.1,1X,I10,1X,I10)') 'Total:',XXBAST,NNBSTRT(1),NNBSTRT(2)
  write(LUPRI,'(A)') '----------------------------------------'
  write(LUPRI,'(/,A,I10,A,I10,A,F7.2,A)') 'First Reduced Set:',NSHP(1),' of',NNSHL_TOT,' shell pairs contribute (',PCT1,'%)'
  write(LUPRI,'(A,I10,A,I10,A,F7.2,A)') 'Curr. Reduced Set:',NSHP(2),' of',NNSHL_TOT,' shell pairs contribute (',PCT2,'%)'
end if

end subroutine CHO_PRTRED
