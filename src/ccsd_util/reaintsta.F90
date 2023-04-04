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

subroutine reaintsta(wrk,wrksize)
! this routine reads integral file INTSTA (reorg), which contains
! following integrals: foka,fokb,
! <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
! <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
! <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
!
! two electron integrals are read to their fix files,
! foka,fokb are read to N,P help files
!
! use and destroy : V1-3, N,P

use ccsd_global, only: daddr, iokey, n, p, v1, v2, v3, w01, w02, w03, w11, w12, w13, w14
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: f_iostat, f_recl, lunsta, rc
logical(kind=iwp) :: is_error

!* open INTSTA file
lunsta = 1
if (iokey == 1) then
  !  Fortran IO
  call molcas_open_ext2(lunsta,'INTSTA','sequential','unformatted',f_iostat,.false.,f_recl,'unknown',is_error)
  !open(unit=lunsta,file='INTSTA',form='unformatted')

else
  !  MOLCAS IO
  call daname(lunsta,'INTSTA')
  daddr(lunsta) = 0
end if

!1 read foka to N
call getmediate(wrk,wrksize,lunsta,n%pos0,n%d,n%i,rc)

!2 read fokb to P
call getmediate(wrk,wrksize,lunsta,p%pos0,p%d,p%i,rc)

!3 read <kl||ij>aaaa to W01
call getmediate(wrk,wrksize,lunsta,w01%pos0,w01%d,w01%i,rc)

!4 read <kl||ij>bbbb to W02
call getmediate(wrk,wrksize,lunsta,w02%pos0,w02%d,w02%i,rc)

!5 read <kl||ij>abab to W03
call getmediate(wrk,wrksize,lunsta,w03%pos0,w03%d,w03%i,rc)

!6 read <ie||mn>aaaa to W11
call getmediate(wrk,wrksize,lunsta,w11%pos0,w11%d,w11%i,rc)

!7 read <ie||mn>bbbb to W12
call getmediate(wrk,wrksize,lunsta,w12%pos0,w12%d,w12%i,rc)

!8 read <ie||mn>abab to W13
call getmediate(wrk,wrksize,lunsta,w13%pos0,w13%d,w13%i,rc)

!9 read <ie||mn>baab to W14
call getmediate(wrk,wrksize,lunsta,w14%pos0,w14%d,w14%i,rc)

!10 read <ab||ij>aaaa to V1
call getmediate(wrk,wrksize,lunsta,v1%pos0,v1%d,v1%i,rc)

!11 read <ab||ij>bbbb to V2
call getmediate(wrk,wrksize,lunsta,v2%pos0,v2%d,v2%i,rc)

!12 read <ab||ij>abab to V3
call getmediate(wrk,wrksize,lunsta,v3%pos0,v3%d,v3%i,rc)

!* close INTSTA file

if (iokey == 1) then
  ! Fortran IO
  close(lunsta)

else
  ! MOLCAS IO
  call daclos(lunsta)
end if

return

end subroutine reaintsta
