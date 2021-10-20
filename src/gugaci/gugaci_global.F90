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

module gugaci_global

use Constants, only: Two, Three, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"

integer(kind=iwp), parameter :: lenintegral = 4, loputmp = 10000, max_atom = 200, max_extorb = 300, max_h0 = 13000, &
                                max_innorb = 100, max_iter = 100, max_kspace = 40, max_lpext_mode = 60000, max_node = 36000, &
                                max_orb = 500, max_ref = 128, max_root = 20, max_tmpvalue = 1000000, max_vector = 800000000, &
                                max_wei = 208000, maxpl = 300000, mtmp = max_orb*(max_orb-1)/2, ntrabuf = 9600, ntratoc = 106, &
                                nvaltype = 200000
integer(kind=iwp), parameter :: istep_occ(4) = [0,1,1,2], &
                                map_jplr(6,6) = reshape([25,23,17,10,24,18, &  !v
                                                         26,19,13, 6,22, 0, &  !d
                                                          0,14,11, 2, 0, 0, &  !t
                                                          0, 7, 3, 1, 9, 5, &  !s
                                                          0,21, 0, 8,20,15, &  !dd
                                                          0, 0, 0, 4,16,12], & !tt
                                                        shape(map_jplr))
real(kind=wp), parameter :: v_onevsqtwo = sqrt(Half), v_sqthree = sqrt(Three), v_sqthreevsqtwo = sqrt(OneHalf), v_sqtwo = sqrt(Two)

integer(kind=iwp) :: ibsm_ext(mxSym), ican_a(max_orb), ican_b(mtmp+max_orb), icano_nnend, icano_nnsta, icnt_base, idisk_array(13), &
                     idisk_lp, idownwei_g131415, iesm_ext(mxSym), ifrno(max_h0), ildownwei_segdd, ilsegdownwei, iml, imr, &
                     index_lpext3(max_extorb,max_innorb), index_lpext4(max_extorb,max_innorb), index_lpext5(max_extorb), &
                     indx(max_kspace), int_dd_drl, int_dd_offset(8,8), ip2_aa_ext_base, ip2_dd_ext_base, ip3_abd_ext_base, &
                     ip4_abcd_ext_base(180), ipae, ipael, ipaety, irdownwei_segdd, iref_occ(max_innorb,max_ref), irf, &
                     irfno(max_ref), irsegdownwei, iseg_downwei(25), iseg_sta(26), iseg_upwei(25), isegdownwei, isegsta, &
                     isegupwei, ism_g1415, ism_g2g4, ivaluesta_g26, iw_downwei(41,25), iw_sta(41,25), iweista_g25, iweista_g26, &
                     iweista_g28, iwt_orb_ext(max_extorb,max_extorb), iwt_sm_s_ext, jb_sys, jd(8), jml, jmr, jp2(max_orb), &
                     jp3(max_orb), jpad, jpad_upwei(82), jpadl, jpadlr, jpae, jpae_downwei(25), jpael, jpel, jper, jph_, &
                     jroute_sys, js(8), jt(8), jud(max_innorb), just(max_innorb,max_innorb), jv, jwl, jwr, lenvec, line, linelp, &
                     log_prod, lp_count(22), lpblock, lpblock_dd, lpblock_ds, lpblock_dt, lpblock_dv, lpblock_sd, lpblock_ss, &
                     lpblock_st, lpblock_sv, lpblock_td, lpblock_ts, lpblock_tt, lpblock_tv, lpblock_vd, lpend34a, lpend34b, &
                     lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei(max_lpext_mode*3), lpsta34a, lpsta34b, lpsta35a, lpsta35b, &
                     lpsta36a, lpsta36b, lrg, lrs, lsm(max_orb), lsm_inn(max_innorb), lsmorb(max_orb), LuCiDen, LuCiDia, LuCiInt, &
                     LuCiMO, LuCiTv1, LuCiTv2, LuCiVec, LuDrt, LuLoop, LuOneMO, LuTwoMO, m_jc, m_jd, map_orb_order(max_orb), &
                     maxciiter, maxintseg, mcroot, mhlp, mhlpmax, mhsum, mjn(2*max_root), mroot, mth_eigen, mtype, mxnode, n_ref, &
                     nabc, nci_dim, nci_h0, ncibl(max_orb), ncibl_all(max_orb), ndim, ndim_h0, ndr, ng_sm, ngw2(max_orb), &
                     ngw3(max_orb), ngw4(max_orb), nint_g25, nint_g28, nlg1, nlg2, nlsm_all(mxSym), nlsm_bas(mxSym), &
                     nlsm_dbl(mxSym), nlsm_ext(mxSym), nlsm_frz(mxSym), no(0:max_innorb), nohy, noidx(8), norb_act, norb_all, &
                     norb_dbl, norb_dz, norb_ext, norb_frz, norb_inn, norb_number(max_orb), np3_abd_ext, ns_sm, nstart_act, &
                     nstaval(nvaltype), nu_ad(41), nu_ae(25), nvalue(nvaltype), nvalue_space_ss, nwalk(0:max_orb), nwei_g25, &
                     nwei_g26, nwei_g28 !, jpadlrel(26), len_str, naorbs, ndjgrop, ndjmod, numat
real(kind=wp) :: cm_cri, dm1tmp(max_orb*(max_orb+1)/2), ecih0(max_root), escf(max_root), fg, pd, pdd, pror, ps1, ps2, ps3, ps4, &
                 pt, ptt, spin, value_lpext3(max_extorb,max_innorb), value_lpext4(max_extorb,max_innorb), &
                 value_lpext5(max_extorb), vd(max_kspace), vdint(0:max_orb,0:max_orb), ve(max_kspace), &
                 viasum_0(max_innorb,max_extorb), viasum_1(max_innorb,max_extorb), vijkk_0sum(max_tmpvalue), &
                 vijkk_1sum(max_tmpvalue), voint(0:max_orb,0:max_orb), vp(max_kspace*(max_kspace+1)/2), vpotnuc, vthrealp, &
                 vthreen, vthreresid, vu(max_kspace,max_kspace), w0, w0_d1d(2), w0_d1d1(3), w0_d1s(4), w0_d1t1, w0_d1v(2), &
                 w0_dd(3), w0_dd1, w0_ds(3), w0_dt, w0_dv(2), w0_plp, w0_sd(16), w0_sd1(13), w0_sdplp, w0_sdplp25, w0_ss(20), &
                 w0_sv(3), w0_t1d1(5), w0_t1t1(3), w0_td(5), w0_tt(3), w0_vv, w0g13a, w0g14a, w0g15a, w0g25, w0g25a, w0g25b, &
                 w0g26a, w0g26b, w0g27, w0g28a, w0g28b, w0g29, w0g2a, w0g2b, w0g30, w0g31, w0g32, w0g34a, w0g34b, w0g35a, w0g35b, &
                 w0g36a, w0g36b, w0g4a, w0g4b, w0gdd, w0plp25, w0plp26, w0plp27, w0plp28, w0plp29, w0plp30, w0plp31, w0plp32, w1, &
                 w1_d1d(2), w1_d1d1(3), w1_d1s(4), w1_d1t1, w1_d1v(2), w1_dd(3), w1_dd1, w1_ds(3), w1_dt, w1_plp, w1_sd(16), &
                 w1_sd1(13), w1_sdplp, w1_sdplp25, w1_ss(20), w1_st(7), w1_st1(4), w1_sv(3), w1_t1d1(5), w1_t1s(7), w1_t1t1(3), &
                 w1_t1v, w1_td(5), w1_ts(4), w1_tt(3), w1_tv, w1g14a, w1g15a, w1g25a, w1g25b, w1g26a, w1g26b, w1g27, w1g28a, &
                 w1g28b, w1g2a, w1g2b, w1g31, w1g32, w1g34a, w1g34b, w1g35a, w1g35b, w1g36a, w1g36b, w1g4a, w1g4b, w1gdd, w1plp27, &
                 w1plp31, w1plp32
                 !, cf(max_orb,max_orb), dm1(max_orb,max_orb), dxyz(3,max_atom), p(max_orb,max_orb), xlgrn(max_orb,max_orb)
logical(kind=iwp) :: logic_assign_actorb, logic_calpro, logic_dh, logic_g13, logic_g1415, logic_g25a, logic_g25b, logic_g26, &
                     logic_g28a, logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                     logic_g49a, logic_g49b, logic_g50, logic_grad, logic_inivec_read, logic_mr !, logic_mrelcas , logic_tdav
!character(len=256) :: tmpdir
character(len=8) :: FnOneMO, FnTwoMO
integer(kind=iwp), allocatable :: ihy(:), ihyl(:), index_lpext(:), index_lpext1(:), index_lpext2(:), intind_abkk(:), &
                                  intind_iabc(:), intind_iaqq(:), intind_ijab(:), intind_ijcc(:), intind_ijka(:), &
                                  intspace_abkk(:), intspace_ijab(:), intspace_ijcc(:), iy(:,:), iyl(:,:), ja(:), jb(:), jeh(:), &
                                  jj(:,:), jj_sub(:,:), jjl_sub(:,:), jm(:), jph(:), jphy(:), jphyl(:), jwh(:), kk(:), loij(:), &
                                  loij_all(:), loijk(:), loijk_all(:), lp_coe(:,:), lp_head(:), lp_ltail(:), lp_lwei(:), &
                                  lp_rtail(:), lp_rwei(:), lpnew_coe(:,:), lpnew_head(:), lpnew_ltail(:), lpnew_lwei(:), &
                                  lpnew_rtail(:), lpnew_rwei(:)
real(kind=wp), allocatable :: denm1(:), denm2(:), th(:), thh(:), value_lpext(:), value_lpext1(:), value_lpext2(:), vcm(:), &
                              vector1(:), vector2(:), vint_ci(:), vplp_w0(:), vplp_w1(:), vplpnew_w0(:), vplpnew_w1(:)
logical(kind=iwp), allocatable :: logic_br(:), logic_newbr(:)

public :: cm_cri, denm1, denm2, dm1tmp, ecih0, escf, fg, FnOneMO, FnTwoMO, ibsm_ext, ican_a, ican_b, icano_nnend, icano_nnsta, &
          icnt_base, idisk_array, idisk_lp, idownwei_g131415, iesm_ext, ifrno, ihy, ihyl, ildownwei_segdd, ilsegdownwei, iml, imr, &
          index_lpext, index_lpext1, index_lpext2, index_lpext3, index_lpext4, index_lpext5, indx, int_dd_drl, int_dd_offset, &
          intind_abkk, intind_iabc, intind_iaqq, intind_ijab, intind_ijcc, intind_ijka, intspace_abkk, intspace_ijab, &
          intspace_ijcc, ip2_aa_ext_base, ip2_dd_ext_base, ip3_abd_ext_base, ip4_abcd_ext_base, ipae, ipael, ipaety, &
          irdownwei_segdd, iref_occ, irf, irfno, irsegdownwei, iseg_downwei, iseg_sta, iseg_upwei, isegdownwei, isegsta, &
          isegupwei, ism_g1415, ism_g2g4, istep_occ, ivaluesta_g26, iw_downwei, iw_sta, iweista_g25, iweista_g26, iweista_g28, &
          iwt_orb_ext, iwt_sm_s_ext, iy, iyl, ja, jb, jb_sys, jd, jeh, jj, jj_sub, jjl_sub, jm, jml, jmr, jp2, jp3, jpad, &
          jpad_upwei, jpadl, jpadlr, jpae, jpae_downwei, jpael, jpel, jper, jph, jph_, jphy, jphyl, jroute_sys, js, jt, jud, just, &
          jv, jwh, jwl, jwr, kk, lenintegral, lenvec, line, linelp, log_prod, logic_assign_actorb, logic_br, logic_calpro, &
          logic_dh, logic_g13, logic_g1415, logic_g25a, logic_g25b, logic_g26, logic_g28a, logic_g2g4a, logic_g2g4b, logic_g34a, &
          logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, logic_g49a, logic_g49b, logic_g50, logic_grad, &
          logic_inivec_read, logic_mr, logic_newbr, loij, loij_all, loijk, loijk_all, loputmp, lp_coe, lp_count, lp_head, &
          lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpblock, lpblock_dd, lpblock_ds, lpblock_dt, lpblock_dv, lpblock_sd, lpblock_ss, &
          lpblock_st, lpblock_sv, lpblock_td, lpblock_ts, lpblock_tt, lpblock_tv, lpblock_vd, lpend34a, lpend34b, lpend35a, &
          lpend35b, lpend36a, lpend36b, lpext_wei, lpnew_coe, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, &
          lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, lpsta36b, lrg, lrs, lsm, lsm_inn, lsmorb, LuCiDen, LuCiDia, LuCiInt, &
          LuCiMO, LuCiTv1, LuCiTv2, LuCiVec, LuDrt, LuLoop, LuOneMO, LuTwoMO, m_jc, m_jd, map_jplr, map_orb_order, max_extorb, &
          max_h0, max_innorb, max_iter, max_kspace, max_node, max_orb, max_ref, max_root, max_tmpvalue, max_vector, max_wei, &
          maxciiter, maxintseg, maxpl, mcroot, mhlp, mhlpmax, mhsum, mjn, mroot, mth_eigen, mtype, mxnode, n_ref, nabc, nci_dim, &
          nci_h0, ncibl, ncibl_all, ndim, ndim_h0, ndr, ng_sm, ngw2, ngw3, ngw4, nint_g25, nint_g28, nlg1, nlg2, nlsm_all, &
          nlsm_bas, nlsm_dbl, nlsm_ext, nlsm_frz, no, nohy, noidx, norb_act, norb_all, norb_dbl, norb_dz, norb_ext, norb_frz, &
          norb_inn, norb_number, np3_abd_ext, ns_sm, nstart_act, nstaval, ntrabuf, ntratoc, nu_ad, nu_ae, nvaltype, nvalue, &
          nvalue_space_ss, nwalk, nwei_g25, nwei_g26, nwei_g28, pd, pdd, pror, ps1, ps2, ps3, ps4, pt, ptt, spin, th, thh, &
          v_onevsqtwo, v_sqthree, v_sqthreevsqtwo, v_sqtwo, value_lpext, value_lpext1, value_lpext2, value_lpext3, value_lpext4, &
          value_lpext5, vcm, vd, vdint, ve, vector1, vector2, viasum_0, viasum_1, vijkk_0sum, vijkk_1sum, vint_ci, voint, vp, &
          vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, vpotnuc, vthrealp, vthreen, vthreresid, vu, w0, w0_d1d, w0_d1d1, w0_d1s, &
          w0_d1t1, w0_d1v, w0_dd, w0_dd1, w0_ds, w0_dt, w0_dv, w0_plp, w0_sd, w0_sd1, w0_sdplp, w0_sdplp25, w0_ss, w0_sv, w0_t1d1, &
          w0_t1t1, w0_td, w0_tt, w0_vv, w0g13a, w0g14a, w0g15a, w0g25, w0g25a, w0g25b, w0g26a, w0g26b, w0g27, w0g28a, w0g28b, &
          w0g29, w0g2a, w0g2b, w0g30, w0g31, w0g32, w0g34a, w0g34b, w0g35a, w0g35b, w0g36a, w0g36b, w0g4a, w0g4b, w0gdd, w0plp25, &
          w0plp26, w0plp27, w0plp28, w0plp29, w0plp30, w0plp31, w0plp32, w1, w1_d1d, w1_d1d1, w1_d1s, w1_d1t1, w1_d1v, w1_dd, &
          w1_dd1, w1_ds, w1_dt, w1_plp, w1_sd, w1_sd1, w1_sdplp, w1_sdplp25, w1_ss, w1_st, w1_st1, w1_sv, w1_t1d1, w1_t1s, &
          w1_t1t1, w1_t1v, w1_td, w1_ts, w1_tt, w1_tv, w1g14a, w1g15a, w1g25a, w1g25b, w1g26a, w1g26b, w1g27, w1g28a, w1g28b, &
          w1g2a, w1g2b, w1g31, w1g32, w1g34a, w1g34b, w1g35a, w1g35b, w1g36a, w1g36b, w1g4a, w1g4b, w1gdd, w1plp27, w1plp31, w1plp32

end module gugaci_global
