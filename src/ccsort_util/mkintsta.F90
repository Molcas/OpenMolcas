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

subroutine mkintsta(wrk,wrksize,foka,fokb)
! this routine produces integral file INTSTA, which contains
! following integrals: foka,fokb,
! <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
! <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
! <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
!
! N.B. 1. work file #1 is used for <ij|pq> integrals, #2,3,4
! must be free. posb0 must be defined
! N.B. 2. this routine can be used only after definition of <ij|pq>
! N.B. 3. this routine use following help routines:
! expandfok
! dawrtmediate (from SYMM)

use ccsort_global, only: daddr, iokey, map2
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: foka(*), fokb(*)
integer(kind=iwp) :: lunsta, rc

! open INTSTA file
lunsta = 21
if (iokey == 1) then
  ! Fortarn IO
  call molcas_binaryopen_vanilla(lunsta,'INTSTA')
  !open(unit=lunsta,file='INTSTA',form='unformatted')

else
  ! MOLCAS IO
  call daname(lunsta,'INTSTA')
  daddr(lunsta) = 0
end if

! expand foka into work #2 and write to INTSTA
call expandfok(wrk,wrksize,foka)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! expand fokb into work #2 and write to INTSTA
call expandfok(wrk,wrksize,fokb)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <kl||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,4,1,1,1,1,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <kl||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,4,2,2,2,2,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <kl| ij>abab from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,0,1,2,1,2,1,0)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ka||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,3,1,3,1,1,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ka||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,3,2,4,2,2,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ka| ij>abab from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,0,1,4,1,2,1,0)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ka||ij>baab from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,0,2,3,1,2,0,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ab||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,4,3,3,1,1,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ab||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,4,4,4,2,2,1,1)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! get #2 <ab| ij>abab from #1 <p,q|i,j> and write to INTSTA
call exppqij(wrk,wrksize,0,3,4,1,2,1,0)
call dawrtmediate(wrk,wrksize,lunsta,map2,rc)

! close INTSTA file

if (iokey == 1) then
  ! Fortran IO
  close(lunsta)

else
  ! MOLCAS IO
  call daclos(lunsta)
end if

return

end subroutine mkintsta
