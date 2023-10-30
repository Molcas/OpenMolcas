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

subroutine hess_svb1_cvb(orbs,civecp,civbs,civb,citmp,orbinv,sorbs,owrk,dvbdet,grad1,grad2,hessorb,vec1,iorts,hessinp,hessout)

use casvb_global, only: aa1, aa2, n_cihess, n_orbhess, ndet, ndetvb, nfrag, norb, nort, npr, nprorb, nprvb, nvb, oaa2, oaa3, proj, &
                        projcas, strucopt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: orbs(norb,norb), civb(0:ndet), orbinv(norb,norb), grad1(npr), grad2(npr), hessorb(nprorb,nprorb)
real(kind=wp), intent(inout) :: civecp(0:ndet), civbs(0:ndet), citmp(0:ndet), dvbdet(ndetvb), hessinp(npr)
real(kind=wp), intent(out) :: sorbs(norb,norb), owrk(norb,norb), hessout(npr)
! VEC1 dimension is MAX(NPRORB,NDETVB)
real(kind=wp), intent(_OUT_) :: vec1(*)
integer(kind=iwp), intent(in) :: iorts(2,nort)
integer(kind=iwp) :: iorb, iort, jorb, ki, kj, korb, lj, lorb
real(kind=wp) :: corr1, fac1, fac2, g1f, g2f, hess_ci_nrm, hess_orb_nrm
logical(kind=iwp) :: orbopt2, strucopt2
real(kind=wp), allocatable :: cvb(:), cvbdet(:), owrk2(:,:), owrk3(:,:)
real(kind=wp), external :: ddot_, dnrm2_

hess_orb_nrm = dnrm2_(nprorb,hessinp,1)
orbopt2 = hess_orb_nrm > 1.0e-10_wp
if (nprvb > 0) then
  hess_ci_nrm = dnrm2_(nprvb,hessinp(nprorb+1:),1)
else
  hess_ci_nrm = Zero
end if
strucopt2 = strucopt .and. (hess_ci_nrm > 1.0e-10_wp)
if (orbopt2 .and. (.not. strucopt2)) n_orbhess = n_orbhess+1
if (strucopt2 .and. (.not. orbopt2)) n_cihess = n_cihess+1

call trnsps(norb,norb,orbs,owrk)
call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)

hessout(:) = Zero
if (orbopt2) call mxatb_cvb(hessorb,hessinp,nprorb,nprorb,1,hessout)
! Combinations of gradients:
g1f = ddot_(npr,grad1,1,hessinp,1)
g2f = ddot_(npr,grad2,1,hessinp,1)
fac1 = g1f*oaa3+g2f*aa2
fac2 = g1f*aa2
hessout(:) = hessout(:)+fac1*grad1(:)+fac2*grad2(:)

if (orbopt2 .and. strucopt) then
  call mxunfold_cvb(hessinp,owrk,norb)
  call dgetmi(owrk,norb,norb)
  call cizero_cvb(citmp)
  call oneexc_cvb(civecp,citmp,owrk,.true.,2)
  call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.false.)
  hessout(nprorb+1:nprorb+nprvb) = hessout(nprorb+1:nprorb+nprvb)+aa1*vec1(nprorb+1:nprorb+nprvb)
  call mma_allocate(owrk2,norb,norb,label='owrk2')
  call mma_allocate(owrk3,norb,norb,label='owrk3')
  call trnsps(norb,norb,owrk,owrk2)
  if (.not. (proj .or. projcas)) then
    call mxatb_cvb(sorbs,owrk2,norb,norb,norb,owrk3)
    call mxatb_cvb(owrk3,orbinv,norb,norb,norb,owrk2)
    owrk(:,:) = owrk(:,:)+owrk2(:,:)
  end if
  call mma_deallocate(owrk2)
  call mma_deallocate(owrk3)
  call cizero_cvb(citmp)
  call oneexc_cvb(civbs,citmp,owrk,.true.,2)
  call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.false.)
  hessout(nprorb+1:nprorb+nprvb) = hessout(nprorb+1:nprorb+nprvb)+oaa2*vec1(nprorb+1:nprorb+nprvb)
end if
if (strucopt2) then
  call str2vbc_cvb(hessinp(1+nprorb:),dvbdet)
  call vb2cif_cvb(dvbdet,citmp)
  ! Structure coeff. <-> orbital
  call mkgrd_cvb(citmp,civbs,vec1,dvbdet,nprorb,.true.)
  hessout(1:nprorb) = hessout(1:nprorb)+oaa2*vec1(1:nprorb)
  call mkgrd_cvb(citmp,civecp,vec1,dvbdet,nprorb,.true.)
  hessout(1:nprorb) = hessout(1:nprorb)+aa1*vec1(1:nprorb)
  if (proj .or. projcas) call oneexc_cvb(civb,citmp,hessinp,.false.,1)
  ! Structure coeff. <-> all
  call applyts_cvb(citmp,orbs)
  call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.true.)
  hessout(:) = hessout(:)+oaa2*vec1(1:npr)
  ! 2nd-order term for structure coefficients
  if (nfrag > 1) then
    call str2vbc_cvb(hessinp(1+nprorb),dvbdet)
    call mma_allocate(cvbdet,ndetvb,label='cvbdet')
    call mma_allocate(cvb,nvb,label='cvb')
    call ci2ordr_cvb(civbs,dvbdet,cvbdet)
    call vb2strg_cvb(cvbdet,cvb)
    hessout(nprorb+1:nprorb+nvb) = hessout(nprorb+1:nprorb+nvb)+oaa2*cvb(:)
    call ci2ordr_cvb(civecp,dvbdet,cvbdet)
    call vb2strg_cvb(cvbdet,cvb)
    hessout(nprorb+1:nprorb+nvb) = hessout(nprorb+1:nprorb+nvb)+aa1*cvb(:)
    call mma_deallocate(cvb)
    call mma_deallocate(cvbdet)
  end if
else if (proj .or. projcas) then
  call cizero_cvb(citmp)
  call oneexc_cvb(civb,citmp,hessinp,.false.,1)
  ! Structure coeff. <-> all
  call applyts_cvb(citmp,orbs)
  call mkgrd_cvb(civb,citmp,vec1,dvbdet,npr,.true.)
  hessout(:) = hessout(:)+oaa2*vec1(1:npr)
end if

if (orbopt2 .and. (nort > 0)) then
  ! Non-linear correction for orthogonality constraints:
  owrk(:,:) = sorbs(:,:)
  call mxinv_cvb(owrk,norb)
  do iort=1,nort
    iorb = iorts(1,iort)
    jorb = iorts(2,iort)
    corr1 = Zero
    do korb=1,norb
      ki = korb+(iorb-1)*(norb-1)
      if (korb > iorb) ki = ki-1
      kj = korb+(jorb-1)*(norb-1)
      if (korb > jorb) kj = kj-1
      if (korb /= iorb) corr1 = corr1+owrk(jorb,korb)*(aa1*grad2(ki)+Half*oaa2*grad1(ki))
      if (korb /= jorb) corr1 = corr1+owrk(iorb,korb)*(aa1*grad2(kj)+Half*oaa2*grad1(kj))
    end do
    corr1 = -Half*corr1
    do korb=1,norb
      if (korb == iorb) cycle
      ki = korb+(iorb-1)*(norb-1)
      if (korb > iorb) ki = ki-1
      do lorb=1,norb
        if (lorb == jorb) cycle
        lj = lorb+(jorb-1)*(norb-1)
        if (lorb > jorb) lj = lj-1
        hessout(ki) = hessout(ki)+sorbs(korb,lorb)*corr1*hessinp(lj)
        hessout(lj) = hessout(lj)+sorbs(korb,lorb)*corr1*hessinp(ki)
      end do
    end do
  end do
end if

return

end subroutine hess_svb1_cvb
