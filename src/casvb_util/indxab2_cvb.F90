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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine indxab2_cvb(indxa,indxb,nstra,nstrb,iocc,nsa,nsb)

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension indxa(nsa), indxb(nsb)
dimension nstra(mxirrep), nstrb(mxirrep)
dimension iocc(norb+1)

call izero(nstra,mxirrep)
call izero(nstrb,mxirrep)
inda = 0
indb = 0
do iisym=1,mxirrep

  call loopstr0_cvb(iocc,index,nalf,norb)
  do
    irp = 1
    do ia=1,nalf
      irp = md2h(irp,ityp(iocc(ia)))
    end do
    if (irp == iisym) then
      inda = inda+1
      nstra(iisym) = nstra(iisym)+1
      indxa(inda) = index
    end if
    call loopstr_cvb(iocc,index,nalf,norb)
    if (index == 1) exit
  end do

  call loopstr0_cvb(iocc,index,nbet,norb)
  do
    irp = 1
    do ib=1,nbet
      irp = md2h(irp,ityp(iocc(ib)))
    end do
    if (irp == iisym) then
      indb = indb+1
      nstrb(iisym) = nstrb(iisym)+1
      indxb(indb) = index
    end if
    call loopstr_cvb(iocc,index,nbet,norb)
    if (index == 1) exit
  end do

end do

return

end subroutine indxab2_cvb
