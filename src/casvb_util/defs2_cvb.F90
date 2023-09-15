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

subroutine defs2_cvb(ifxorb)

implicit real*8(a-h,o-z)
parameter(iunset=-1357924680,unset=-1357924680d0)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "frag_cvb.fh"
dimension ifxorb(mxorb_cvb)

if (icrit == iunset) then
  if ((.not. variat) .and. (imethod /= 6)) then
    icrit = 1
  else
    icrit = 2
  end if
end if
! Set ifxorb
nfxorb = 0
do i=1,norb
  if (ifxorb(i) == 1) nfxorb = nfxorb+1
end do
call izero(ifxorb(norb+1),mxorb_cvb-norb)

! Set STRUCOPT:
if (projcas) then
  strucopt = .false.
else if (imethod == 11) then
  strucopt = .false.
else if (nvb == 1) then
  strucopt = .false.
else
  nfxvbr = nfxvb
  if (lfxvb == 1) nfxvbr = nvb-nfxvb
  nzrvbr = nzrvb
  if (lzrvb == 1) nzrvbr = nvb-nzrvb
  if (nfxvbr >= nvb) then
    strucopt = .false.
  else if (nfxvbr+nzrvbr >= nvb) then
    write(6,*) ' Should check!'
    call abend_cvb()
  else
    strucopt = .true.
  end if
end if
if (isaddle == iunset) isaddle = 0
! If unset, set default for IMETHOD:
if (imethod == iunset) then
  if (isaddle == 0) then
    imethod = 10
  else
    imethod = 7
  end if
  if ((nfxorb == norb) .and. strucopt .and. (nfrag <= 1)) imethod = 4
else if ((isaddle /= 0) .and. (imethod == 1)) then
  call abend_cvb()
end if
! If unset, set default for MXITER:
if (mxiter == iunset) then
  if (imethod /= 4) then
    mxiter = 50
  else
    mxiter = 200
  end if
end if
if (ishstruc == iunset) ishstruc = 0
if (npcf == iunset) npcf = -2
if (ivbweights == iunset) ivbweights = -1
if (iciweights == iunset) iciweights = 0
call tunedefs2_cvb(imethod,.false.)

return

end subroutine defs2_cvb
