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

subroutine makecivb_cvb(civec,civb,cvbdet,orbs,cvb,ic)
! Construct CIVB and CVBDET:
! IC=0 : CIVB will contain full set of structures (if PROJCAS).
! IC=1 : CIVB will contain only VB structures.

use casvb_global, only: gjorb_type
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec(0:ndet), civb(0:ndet), cvbdet(ndetvb), orbs(norb,norb), cvb(nvb)
integer(kind=iwp) :: ic
type(gjorb_type) :: gjorb
real(kind=wp), allocatable :: orbinv(:,:)
logical(kind=iwp), external :: tstcnt_cvb ! ... Content of CI vectors ...

if (tstcnt_cvb(civb,3-ic)) return

if (.not. projcas) then
  call str2vbc_cvb(cvb,cvbdet)
  call vb2cic_cvb(cvbdet,civb)
else if (projcas) then
  call mma_allocate(orbinv,norb,norb,label='orbinv')
  call mma_allocate(gjorb%r,norb,norb,label='gjorb%r')
  call mma_allocate(gjorb%i1,norb,label='gjorb%i1')
  call mma_allocate(gjorb%i2,2,norb*norb,label='gjorb%i2')

  if (memplenty) then
    call getci_cvb(civec)
    call cicopy_cvb(civec,civb)
  else
    call cird_cvb(civb,61001.2_wp)
  end if
  call fmove_cvb(orbs,orbinv,norb*norb)
  call mxinv_cvb(orbinv,norb)
  call gaussj_cvb(orbinv,gjorb)
  call applyt_cvb(civb,gjorb)
  call ci2vbc_cvb(civb,cvbdet)
  call vb2strc_cvb(cvbdet,cvb)
  if (ic == 1) call vb2cic_cvb(cvbdet,civb)

  call mma_deallocate(orbinv)
  call mma_deallocate(gjorb%r)
  call mma_deallocate(gjorb%i1)
  call mma_deallocate(gjorb%i2)
end if
call setcnt_cvb(civb,3-ic)

return

end subroutine makecivb_cvb
