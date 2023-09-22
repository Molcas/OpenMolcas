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

subroutine evb2cas_cvb(orbs,cvb,fx,ioptc,iter)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), fx
integer(kind=iwp) :: ioptc, iter
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, idum
real(kind=wp) :: dx_amx, dxnrm
integer(kind=iwp), external :: mstackr_cvb
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: tstfile_cvb ! ... Files/Hamiltonian available ...

if (tstfile_cvb(66000.2_wp)) then
  i1 = mstackr_cvb(norb*norb+nvb)
  call rdr_cvb(work(i1),norb*norb+nvb,66000.2_wp,0)
  call subvec(work(i1),orbs,work(i1),norb*norb)
  call subvec(work(norb*norb+i1),cvb,work(norb*norb+i1),nvb)
  dxnrm = dnrm2_(norb*norb+nvb,work(i1),1)
  call findamx_cvb(work(i1),norb*norb+nvb,dx_amx,idum)
  call mfreer_cvb(i1)
end if
call wrr_cvb(orbs,norb*norb,66000.2_wp,0)
call wrr_cvb(cvb,nvb,66000.2_wp,norb*norb)
call evb2cas2_cvb(orbs,cvb,ioptc,iter,fx,dxnrm,dx_amx,work(lc(1)),work(lc(2)),work(lc(3)),work(lc(4)),work(lc(5)),work(lw(2)), &
                  work(lw(3)))

return

end subroutine evb2cas_cvb
