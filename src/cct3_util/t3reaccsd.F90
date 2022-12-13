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

subroutine t3reaccsd(wrk,wrksize,eccsd)
! this routine reads CCSD results, it T1 and T2 amplitudes
! and CCSD energy from saverst file
!
! eccsd - Converged CCSD energy (O)

use CCT3_global, only: daddr, filerst, iokey, t11, t12, t21, t22, t23
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(out) :: eccsd
integer(kind=iwp) :: istatus, lunrst, rc1
real(kind=wp) :: dum(1)

!1 open file savename
lunrst = 1

if (iokey == 1) then
  ! Fortran IO
  !open(unit=lunrst,file=filerst,form='unformatted')
  call molcas_binaryopen_vanilla(lunrst,filerst)

else
  ! MOLCAS IO
  call daname(lunrst,filerst)
  daddr(lunrst) = 0
end if

!2 get T1aa
call cct3_getmediate(wrk,wrksize,lunrst,t11,rc1)

!3 get T1bb
call cct3_getmediate(wrk,wrksize,lunrst,t12,rc1)

!4 get T2aaaa
call cct3_getmediate(wrk,wrksize,lunrst,t21,rc1)

!5 get T2bbbb
call cct3_getmediate(wrk,wrksize,lunrst,t22,rc1)

!6 get T2abab
call cct3_getmediate(wrk,wrksize,lunrst,t23,rc1)

!7 get energy,niter
if (iokey == 1) then
  ! Fortran IO
  read(lunrst,iostat=istatus) eccsd,rc1
  if (istatus < 0) then
    write(u6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
    write(u6,*) ' USE CCSD ENERGY FROM CCSD OUTPUT FILE'
    eccsd = Zero
  end if
else

  ! MOLCAS IO
  call ddafile(lunrst,2,dum,1,daddr(lunrst))
  eccsd = dum(1)
end if

if (iokey == 1) then
  ! Fortran IO
  close(lunrst)

else
  ! MOLCAS IO
  call daclos(lunrst)
end if

return

end subroutine t3reaccsd
