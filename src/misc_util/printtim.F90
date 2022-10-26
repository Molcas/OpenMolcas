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

subroutine printtim(nfld1,nfld2)

use Para_Info, only: nProcs

implicit real*8(a-h,o-z)
#include "timtra.fh"
#include "WrkSpc.fh"

iout = 6
iPL = iPrintLevel(-1)
if (iPL <= 2) return

write(iout,*)
call CollapseOutput(1,'Timing information')
write(iout,*)
if (nfld1 > 0) then
  m = nfld_tim*2
  l = nprocs*m
  call gadgop(work(igatim),l,'+')
  write(iout,1)
  do iproc=0,nprocs-1,10
    n = min(10,nprocs-iproc)
    write(iout,5) (k,k=iproc+1,iproc+n)
    ioff = (igatim-1)+iproc*m
    do j=1,min(nfld_tim,nfld1)
      if (namfld(j) /= ' ') then
        write(iout,10) namfld(j),(work(ioff+j+k*m),k=0,n-1)
      end if
    end do
    ioff = ioff+n*m
  end do
  write(iout,11)
  do iproc=0,nprocs-1,10
    n = min(10,nprocs-iproc)
    write(iout,5) (k,k=iproc+1,iproc+n)
    ioff = (igatim-1)+iproc*m+nfld_tim
    do j=1,min(nfld_tim,nfld1)
      if (namfld(j) /= ' ') then
        write(iout,10) namfld(j),(work(ioff+j+k*m),k=0,n-1)
      end if
    end do
    ioff = ioff+n*m
  end do
  call fzero(work(igatim),l)
end if
if (nfld2 > 0) then
  m = nfld_stat
  l = nprocs*m
  call gadgop(work(igastat),l,'+')
  write(iout,12)
  do iproc=0,nprocs-1,10
    n = min(10,nprocs-iproc)
    write(iout,5) (k,k=iproc+1,iproc+n)
    ioff = (igastat-1)+iproc*m
    do j=1,min(nfld2,nfld_stat)
      if (namstat(j) /= ' ') then
        write(iout,15) namstat(j),(work(ioff+j+k*m),k=0,n-1)
      end if
    end do
    ioff = ioff+n*m
  end do
  call fzero(work(igastat),l)
end if
call CollapseOutput(0,'Timing information')
write(iout,*)

return

1 format(/' CPU times:')
5 format(t20,10i10)
10 format(1x,a20,t21,10f10.2)
11 format(/' Elapsed times:')
12 format(/' Task statistic:')
15 format(1x,a20,t22,10f10.0)

end subroutine printtim
