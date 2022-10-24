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

subroutine t3reaintsta(wrk,wrksize)
! this routine reads integral file INTSTA (reorg), which contains
! following integrals: foka,fokb,
! <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab - naplano
! <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
! <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
!
! two electron integrals are read to their fix files,
! foka,fokb are read to N,P help files
!
! use and destroy : N,P

use CCT3_global, only: daddr, iokey, n, px, w11, w12, w13, w14, w21, w22, w23
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp) :: lunsta, rc

!* open INTSTA file
lunsta = 1
if (iokey == 1) then
  ! Fortran IO
  !open(unit=lunsta,file='INTSTA',form='unformatted')
  call molcas_binaryopen_vanilla(lunsta,'INTSTA')

else
  ! MOLCAS IO
  call daname(lunsta,'INTSTA')
  daddr(lunsta) = 0
end if

!1 read foka to N
call cct3_getmediate(wrk,wrksize,lunsta,n,rc)

!2 read fokb to P
call cct3_getmediate(wrk,wrksize,lunsta,px,rc)

!3 read <kl||ij>aaaa to W23 - naplano
call cct3_getmediate(wrk,wrksize,lunsta,w23,rc)

!4 read <kl||ij>bbbb to W23 - naplano
call cct3_getmediate(wrk,wrksize,lunsta,w23,rc)

!5 read <kl||ij>abab to W23 - naplano
call cct3_getmediate(wrk,wrksize,lunsta,w23,rc)

!6 read <ie||mn>aaaa to W11
call cct3_getmediate(wrk,wrksize,lunsta,w11,rc)

!7 read <ie||mn>bbbb to W12
call cct3_getmediate(wrk,wrksize,lunsta,w12,rc)

!8 read <ie||mn>abab to W13
call cct3_getmediate(wrk,wrksize,lunsta,w13,rc)

!9 read <ie||mn>baab to W14
call cct3_getmediate(wrk,wrksize,lunsta,w14,rc)

!10 read <ab||ij>aaaa to W21
call cct3_getmediate(wrk,wrksize,lunsta,w21,rc)

!11 read <ab||ij>bbbb to W22
call cct3_getmediate(wrk,wrksize,lunsta,w22,rc)

!12 read <ab||ij>abab to W23
call cct3_getmediate(wrk,wrksize,lunsta,w23,rc)

!* close INTSTA file

if (iokey == 1) then
  ! Fortran IO
  close(lunsta)

else
  ! MOLCAS IO
  call daclos(lunsta)
end if

return

end subroutine t3reaintsta
