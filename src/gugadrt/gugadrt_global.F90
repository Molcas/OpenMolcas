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

module gugadrt_global

use Definitions, only: wp, iwp

implicit none
private

#ifndef _I8_
integer(kind=iwp), parameter :: max_ref = 64, iintbit = 32, n32int = 4, n16int = 2
#else
integer(kind=iwp), parameter :: max_ref = 128, iintbit = 64, n32int = 2, n16int = 1
#endif

integer(kind=iwp), parameter :: max_node = 36000, max_orb = 260, max_innorb = 60
integer(kind=iwp), parameter :: mul_tab(8,8) = reshape([1,2,3,4,5,6,7,8, &
                                                        2,1,4,3,6,5,8,7, &
                                                        3,4,1,2,7,8,5,6, &
                                                        4,3,2,1,8,7,6,5, &
                                                        5,6,7,8,1,2,3,4, &
                                                        6,5,8,7,2,1,4,3, &
                                                        7,8,5,6,3,4,1,2, &
                                                        8,7,6,5,4,3,2,1], shape(mul_tab))

real(kind=wp) :: spin
integer(kind=iwp) :: iseg_downwei(25), iseg_sta(26), ja_sys, jb_sys, jc_sys, jpad_upwei(41), jroute_sys, lsm_inn(max_innorb), &
                     mxnode, n_electron, n_ref,  nci_dim, ng_sm, nlsm_all(8), nlsm_bas(8), nlsm_ext(8), nlsmddel(8), nlsmedel(8), &
                     norb_act, norb_all, norb_dbl, norb_dz, norb_ext, norb_frz, norb_inn, ns_sm, nstart_act, nu_ad(41), nu_ae(25)
logical(kind=iwp) :: logic_mr, logic_mrelcas
integer(kind=iwp) :: iref_occ(max_innorb,max_ref)
integer(kind=iwp) :: ja(max_node), jb(max_node), jm(0:max_node), jj(4,0:max_node), kk(0:max_node), no(0:max_innorb), jv, jd(8), &
                     jt(8), js(8)
integer(kind=iwp) :: iprint, ludrt

public :: max_ref, iintbit, n32int, n16int
public :: max_innorb, max_node, max_orb, mul_tab
public :: iseg_downwei, iseg_sta, ja_sys, jb_sys, jc_sys, jpad_upwei, jroute_sys, logic_mr, logic_mrelcas, lsm_inn, mxnode, &
          n_electron, n_ref, nci_dim, ng_sm, nlsm_all, nlsm_bas, nlsm_ext, nlsmddel, nlsmedel, norb_act, norb_all, norb_dbl, &
          norb_dz, norb_ext, norb_frz, norb_inn, ns_sm, nstart_act, nu_ad, nu_ae, spin
public :: iref_occ
public :: ja, jb, jm, jj, kk, no, jv, jd, jt, js
public :: iprint, ludrt

end module gugadrt_global
