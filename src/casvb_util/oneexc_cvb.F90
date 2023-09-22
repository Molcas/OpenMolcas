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

subroutine oneexc_cvb(cfrom,cto,vij,diag,iPvb)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cfrom(*), cto(*), vij(*)
logical(kind=iwp) :: diag
integer(kind=iwp) :: iPvb
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: icfrom, icto, idens, ivij2, nvij
integer(kind=iwp), external :: mstackr_cvb

idens = 0
icfrom = nint(cfrom(1))
icto = nint(cto(1))

if (iform_ci(icfrom) /= 0) then
  write(u6,*) ' Unsupported format in ONEEXC/ONEDENS :',iform_ci(icfrom)
  call abend_cvb()
else if (iform_ci(icto) /= 0) then
  write(u6,*) ' Unsupported format in ONEEXC/ONEDENS :',iform_ci(icto)
  call abend_cvb()
end if

call oneexc2_cvb(work(iaddr_ci(icfrom)),work(iaddr_ci(icto)),vij,iwork(ll(1)),iwork(ll(2)),iwork(ll(5)),iwork(ll(6)),work(ll(9)), &
                 work(ll(10)),iwork(ll(11)),iwork(ll(12)),iwork(ll(13)),iwork(ll(14)),npvb,nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas, &
                 sc,absym(3),diag,idens,iPvb)

! If projcas and iPvb=0 we asume the normal density/1-ex. is required:
if (projcas .and. (iPvb /= 0)) then
  if (diag) then
    nvij = norb*norb
  else
    nvij = norb*(norb-1)
  end if
  ivij2 = mstackr_cvb(nvij)
  if (idens == 0) then
    call fmove_cvb(vij,work(ivij2),nvij)
    call dscal_(nvij,-One,work(ivij2),1)
  else
    call fzero(work(ivij2),nvij)
  end if
  call oneexc2_cvb(work(iaddr_ci(icfrom)),work(iaddr_ci(icto)),work(ivij2),iwork(ll(1)),iwork(ll(2)),iwork(ll(5)),iwork(ll(6)), &
                   work(ll(9)),work(ll(10)),iwork(ll(11)),iwork(ll(12)),iwork(ll(13)),iwork(ll(14)),npvb,nda,ndb,n1a,n1b,nam1, &
                   nbm1,norb,projcas,sc,absym(3),diag,idens,3-iPvb)
  if (idens == 1) call daxpy_(nvij,-One,work(ivij2),1,vij,1)
  call mfreer_cvb(ivij2)
end if

return

end subroutine oneexc_cvb
