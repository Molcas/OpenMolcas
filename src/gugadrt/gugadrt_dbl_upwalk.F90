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

subroutine gugadrt_dbl_upwalk()

use gugadrt_global, only: jpad_upwei, jroute_sys, lsm_inn, mxnode, ng_sm, norb_dbl, norb_dz, norb_frz, nu_ad, ns_sm
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iw, lri, lrj, lsmi, lsmid, lsmij, lsmit, lsmj, no_d, no_t, node

if (norb_dbl == 1) then
  ! v(1),d(2-9),s(18-25)             for s=0
  ! v(1),d(2-9),s(18-25),d'(26-33)   for s<>0

  mxnode = 17+ng_sm
  lri = norb_frz+1
  lsmi = lsm_inn(lri)
  lsmid = Mul(lsmi,ns_sm)
  ! for node_v
  nu_ad(1) = 1
  jpad_upwei(1) = 1
  ! for node_d
  nu_ad(1+lsmid) = 1+lsmid
  jpad_upwei(1+lsmid) = 1
  ! for node_s
  nu_ad(17+ns_sm) = 17+ns_sm
  jpad_upwei(17+ns_sm) = 1

  if (jroute_sys == 1) return
  mxnode = 25+ng_sm
  ! for node_d'
  nu_ad(25+lsmid) = 25+lsmid
  jpad_upwei(25+lsmid) = 1
  return
end if
nu_ad = 0
jpad_upwei = 0

nu_ad(1) = 1
jpad_upwei(1) = 1
if (norb_dbl == 0) then
  mxnode = 1
  return
end if
do lri=norb_frz+1,norb_dz
  lsmi = lsm_inn(lri)
  lsmid = Mul(lsmi,ns_sm)
  no_d = lsmid+1
  jpad_upwei(no_d) = jpad_upwei(no_d)+1
  do lrj=lri+1,norb_dz
    lsmj = lsm_inn(lrj)
    lsmij = Mul(lsmi,lsmj)
    lsmit = Mul(lsmij,ns_sm)
    no_t = lsmit+9
    jpad_upwei(no_t) = jpad_upwei(no_t)+1
  end do
end do
! v(1),d(2-9),t(10-17),s(18-25),d'(26-33),t'(34-41)
select case (jroute_sys)
  case (1)
    mxnode = 25 !v,d,t,s
    jpad_upwei(18:25) = jpad_upwei(10:17)
    jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
  case (2)
    mxnode = 25+8
    jpad_upwei(18:25) = jpad_upwei(10:17)+jpad_upwei(10:17)
    jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
    jpad_upwei(26:33) = jpad_upwei(2:9)
  case (3)
    mxnode = 25+8+8
    jpad_upwei(18:25) = jpad_upwei(10:17)+jpad_upwei(10:17)
    jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
    jpad_upwei(26:33) = jpad_upwei(2:9)
    jpad_upwei(34:41) = jpad_upwei(10:17)
end select
do node=2,mxnode
  iw = jpad_upwei(node)
  if (iw == 0) cycle
  nu_ad(node) = node
end do

return

end subroutine gugadrt_dbl_upwalk
