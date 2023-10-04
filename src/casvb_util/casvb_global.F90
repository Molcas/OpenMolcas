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
! Copyright (C) 2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

module casvb_global

! Format statements depending on iprec:
!--------------------------------------
! form2AD form2AF formAD formAF formChk1 formChk2 formChk3 formcvp formE formMXP1 formMXP2 formMXP3 formMXP4 formMXP5 formMXP6
! formroot formSymW formVBWnorm

use Definitions, only: wp, iwp

implicit none
private

type gjorb_type
  real(kind=wp), allocatable :: r(:,:)
  integer(kind=iwp), allocatable :: i1(:)
  integer(kind=iwp), allocatable :: i2(:,:)
end type gjorb_type

integer(kind=iwp), parameter :: lbuf = 512, mxact_mo = 16, mxdep = 200, mxfield = 500, mxfiles = 40, mxfrag = 10, mxI = 20, &
                                mxirrep_mo = 8, mxMs = 20, mxobj = 100, mxopth = 10, mxprm = 100, mxS = 20, mxstep = 200, &
                                nspinb = 7, nstackrep = 50

integer(kind=iwp) :: i_dep_on_j(mxdep), i2s_fr(mxS,mxfrag), iact_mo(mxact_mo), iaddrm(mxfield), ibuf, ibuffer(lbuf), icase6, &
                     icase7, icnt, icode(mxstep), idan(mxfiles), ifield, ifollow, ifsc_fr(mxfrag), iline = 0, ilv(300), inp, &
                     inputmode, ioffs(mxobj+1), iopt2step(0:30), ioptcode(30), ioptim, ioptstep, ip, ipAnchr, ipdd, ipos(mxstep), &
                     ipp, ipp12e, ipp12s, ipp7, iprint, iprm, iroot, is_set = 0, isaddle, isaddledd, istackrep(nstackrep), iter, &
                     iter12e, iter12s, iter7, izbuffer(lbuf), j_dep_on_i(mxdep), joffs(mxobj+1), joptstep, jroot, lenline, &
                     loopstep, loopstepmx, lstprm(mxprm), maxd, mnion_fr(mxfrag), mxdav, mxion_fr(mxfrag), mxit, mxrhs, n_div, &
                     nact_mo, nalf_fr(mxMs,mxfrag), nbas_mo, nbasf_mo(mxirrep_mo), nbasi_mo(mxirrep_mo), nbasisq_mo, &
                     nbassqf_mo(mxirrep_mo), nbassqi_mo(mxirrep_mo), nbet_fr(mxMs,mxfrag), nbuf, ncnt, nconf_fr(mxfrag), &
                     nconfion_fr(0:mxI,mxfrag), nda_fr(mxMs,mxfrag), ndb_fr(mxMs,mxfrag), ndep_ij, ndep_ji, ndetvb_fr(mxfrag), &
                     ndetvb2_fr(mxfrag), nel_fr(mxfrag), nfield, nfold, nfrag, nfrdim, nline = 0, nlold = 0, nMs_fr(mxfrag), nobj, &
                     nopth1(2), nopth2(2), noptim, noptstep, nortiter, nortiterdd, nparm, nroot, nS_fr(mxfrag), nsym_mo, &
                     nvb_fr(mxfrag), nvbr_fr(mxfrag), nvecmx, nvguess, nvrestart, nvrhs, nvtot, nword
real(kind=wp) :: aa1, aa2, alftol, cnrm, cnrmtol, corenrg, cpropt(mxopth), cvbnrm_fr(mxfrag), delopth1(2), delopth2(2), dfx(6), &
                 dfxmin(2), dfxtol, dx(3,6), eigwrngtol, exp12tol, expct, f1, f2, f3, f4, file_id, fxbest, grd(3,6), grdwrngtol, &
                 hh, hhaccfac(5,2), hhkeep, hhmax(2), hhopt(mxopth), hhrejfac(2), hhstart, hhtol(2), oaa2, oaa3, orththr, &
                 orththrdd, ovraa_try, ovrab, ovrab_try, recn, resthr, resthrdd, safety, sgn(6), signtol, singul(3), ww, ww_try, &
                 zzacclim(4,2), zzmax(6), zzmin(6), zzrejmax(2), zzrejmin(2)
logical(kind=iwp) :: endwhenclose, follow, have_solved_it, maxize, mustdeclare, ndres_ok, release(10), scalesmall(2), up2date(mxobj)
character(len=300) :: line
character(len=20) :: form2AD, form2AF, formAD, formAF, formChk1, formChk2, formChk3, formcvp, formE, formMXP1, formMXP2, formMXP3, &
                     formMXP4, formMXP5, formMXP6, formroot, formSymW, formVBWnorm
character(len=8) :: charobj(mxobj)
type(gjorb_type) :: gjorb, gjorb2, gjorb3
integer(kind=iwp), allocatable :: confsinp(:,:), ia12ind(:), iapr(:), iapr1(:), ib12ind(:), ibpr(:), ibpr1(:), iconfs(:,:), &
                                  idelstr(:), idetvb(:), ifxorb(:), ifxstr(:), ikcoff(:,:,:), iorbrel(:), iorts(:,:), &
                                  ipermzeta(:,:), irels(:,:), irots(:,:), ixapr(:), ixapr1(:), ixbpr(:), ixbpr1(:), izeta(:), &
                                  ndetvbs(:,:), north(:)
integer(kind=iwp), allocatable, target :: i1alf(:,:), i1c(:,:), iafrm(:,:), iato(:,:), icfrm(:,:), icto(:,:), ifnss1(:,:), &
                                          ifnss2(:,:)
real(kind=wp), allocatable :: ap(:,:), axc(:,:), c(:,:), corth(:,:), cvb(:), cvbdet(:), cvbsspn(:), cvbstot(:), cvbtry(:), &
                              dvbdet(:), eigval(:), eigvec(:,:), evbdet(:), grad1(:), grad2(:), gradx(:,:), gsinp(:), &
                              hessorb(:,:), hesst(:,:), odx(:), odxp(:), ograd(:), ogradp(:), orbinv(:,:), orbs(:,:), &
                              orbstry(:,:), owrk(:), owrk2(:,:), relorb(:,:,:), res(:), rhs(:,:), rhsp(:), solp(:), solp_res(:), &
                              sorbs(:,:), span(:,:), sstruc(:,:), sstruc2(:,:), sxc(:,:), symelm(:,:,:), tconstr(:,:), trprm(:,:), &
                              vbdet(:), vec1(:), wdx(:)
real(kind=wp), allocatable, target :: aikcof(:), cikcof(:), civbvecs(:,:), phato(:,:), phcto(:,:)
integer(kind=iwp), pointer :: i1bet(:,:) => null(), ibfrm(:,:) => null(), ibto(:,:) => null()
real(kind=wp), pointer :: bikcof(:) => null(), civb1(:) => null(), civb2(:) => null(), civb3(:) => null(), civb4(:) => null(), &
                          civb5(:) => null(), civb6(:) => null(), civb7(:) => null(), civb8(:) => null(), civbvec(:,:) => null(), &
                          phbto(:,:) => null()

integer(kind=iwp), parameter :: iunset = -1357924680
logical(kind=iwp), parameter :: ifhamil = .true., ifmos = .true.
character(len=*), parameter :: spinb(nspinb) = ['Kotani      ','Serber      ','Rumer       ','Rumer (LT)  ','Projected   ', &
                                                'Determinants','Determinants'], &
                               spinbkw(nspinb) = ['KOTANI  ','SERBER  ','RUMER   ','LTRUMER ','PROJECT ','DET     ', &
                                                  'DETERM  ']

public :: aa1, aa2, aikcof, alftol, ap, axc, bikcof, c, casvb_free, charobj, cikcof,  civb1, civb2, civb3, civb4, civb5, civb6, &
          civb7, civb8, civbvec, civbvecs, cnrm, cnrmtol, confsinp, corenrg, corth, cpropt, cvb, cvbdet, cvbnrm_fr, cvbsspn, &
          cvbstot, cvbtry, delopth1, delopth2, dfx, dfxmin, dfxtol, dvbdet, dx, eigval, eigvec, eigwrngtol, endwhenclose, evbdet, &
          exp12tol, expct, f1, f2, f3, f4, file_id, follow, form2AD, form2AF, formAD, formAF, formChk1, formChk2, formChk3, &
          formcvp, formE, formMXP1, formMXP2, formMXP3, formMXP4, formMXP5, formMXP6, formroot, formSymW, formVBWnorm, fxbest, &
          gjorb, gjorb_type, gjorb2, gjorb3, grad1, grad2, gradx, grd, grdwrngtol, gsinp, have_solved_it, hessorb, hesst, hh, &
          hhaccfac, hhkeep, hhmax, hhopt, hhrejfac, hhstart, hhtol, i_dep_on_j, i1alf, i1bet, i1c, i2s_fr, ia12ind, iact_mo, &
          iaddrm, iafrm, iapr, iapr1, iato, ib12ind, ibfrm, ibpr, ibpr1, ibto, ibuf, ibuffer, icase6, icase7, icfrm, icnt, icode, &
          iconfs, icto, idan, idelstr, idetvb, ifhamil, ifield, ifmos, ifnss1, ifnss2, ifollow, ifsc_fr, ifxorb, ifxstr, ikcoff, &
          iline, ilv, inp, inputmode, ioffs, iopt2step, ioptcode, ioptim, ioptstep, iorbrel, iorts, ip, ipAnchr, ipdd, ipermzeta, &
          ipos, ipp, ipp12e, ipp12s, ipp7, iprint, iprm, irels, iroot, irots, is_set, isaddle, isaddledd, istackrep, iter, &
          iter12e, iter12s, iter7, iunset, ixapr, ixapr1, ixbpr, ixbpr1, izbuffer, izeta, j_dep_on_i, joffs, joptstep, jroot, &
          lbuf, lenline, line, loopstep, loopstepmx, lstprm, maxd, maxize, mnion_fr, mustdeclare, mxact_mo, mxdav, mxdep, mxfield, &
          mxfiles, mxfrag, mxI, mxion_fr, mxirrep_mo, mxit, mxMs, mxobj, mxopth, mxprm, mxrhs, mxS, mxstep, n_div, nact_mo, &
          nalf_fr, nbas_mo, nbasf_mo, nbasi_mo, nbasisq_mo, nbassqf_mo, nbassqi_mo, nbet_fr, nbuf, ncnt, nconf_fr, nconfion_fr, &
          nda_fr, ndb_fr, ndep_ij, ndep_ji, ndetvb_fr, ndetvb2_fr, ndetvbs, ndres_ok, nel_fr, nfield, nfold, nfrag, nfrdim, nline, &
          nlold, nMs_fr, nobj, nopth1, nopth2, noptim, noptstep, north, nortiter, nortiterdd, nparm, nroot, nS_fr, nspinb, &
          nstackrep, nsym_mo, nvb_fr, nvbr_fr, nvecmx, nvguess, nvrestart, nvrhs, nvtot, nword, oaa2, oaa3, odx, odxp, ograd, &
          ogradp, orbinv, orbs, orbstry, orththr, orththrdd, ovraa_try, ovrab, ovrab_try, owrk, owrk2, phato, phbto, phcto, recn, &
          release, relorb, res, resthr, resthrdd, rhs, rhsp, safety, scalesmall, sgn, signtol, singul, solp, solp_res, sorbs, &
          span, spinb, spinbkw, sstruc, sstruc2, sxc, symelm, tconstr, trprm, up2date, vbdet, vec1, wdx, ww, ww_try, zzacclim, &
          zzmax, zzmin, zzrejmax, zzrejmin

contains

subroutine casvb_free()
  use stdalloc, only: mma_deallocate
  if (allocated(civbvecs)) call mma_deallocate(civbvecs)
  nullify(civbvec)
  nullify(civb1)
  nullify(civb2)
  nullify(civb3)
  nullify(civb4)
  nullify(civb5)
  nullify(civb6)
  nullify(civb7)
  nullify(civb8)
  if (allocated(orbinv)) call mma_deallocate(orbinv)
  if (allocated(sorbs)) call mma_deallocate(sorbs)
  if (allocated(owrk2)) call mma_deallocate(owrk2)
  if (allocated(gjorb%r)) call mma_deallocate(gjorb%r)
  if (allocated(gjorb%i1)) call mma_deallocate(gjorb%i1)
  if (allocated(gjorb%i2)) call mma_deallocate(gjorb%i2)
  if (allocated(gjorb2%r)) call mma_deallocate(gjorb2%r)
  if (allocated(gjorb2%i1)) call mma_deallocate(gjorb2%i1)
  if (allocated(gjorb2%i2)) call mma_deallocate(gjorb2%i2)
  if (allocated(gjorb3%r)) call mma_deallocate(gjorb3%r)
  if (allocated(gjorb3%i1)) call mma_deallocate(gjorb3%i1)
  if (allocated(gjorb3%i2)) call mma_deallocate(gjorb3%i2)
  if (allocated(cvbstot)) call mma_deallocate(cvbstot)
  if (allocated(cvbsspn)) call mma_deallocate(cvbsspn)
  if (allocated(cvbdet)) call mma_deallocate(cvbdet)
  if (allocated(dvbdet)) call mma_deallocate(dvbdet)
  if (allocated(evbdet)) call mma_deallocate(evbdet)
  if (allocated(orbstry)) call mma_deallocate(orbstry)
  if (allocated(cvbtry)) call mma_deallocate(cvbtry)
  if (allocated(i1alf)) call mma_deallocate(i1alf)
  nullify(i1bet)
  if (allocated(i1c)) call mma_deallocate(i1c)
  if (allocated(iafrm)) call mma_deallocate(iafrm)
  nullify(ibfrm)
  if (allocated(icfrm)) call mma_deallocate(icfrm)
  if (allocated(iato)) call mma_deallocate(iato)
  nullify(ibto)
  if (allocated(icto)) call mma_deallocate(icto)
  if (allocated(phato)) call mma_deallocate(phato)
  nullify(phbto)
  if (allocated(phcto)) call mma_deallocate(phcto)
  if (allocated(iapr)) call mma_deallocate(iapr)
  if (allocated(ixapr)) call mma_deallocate(ixapr)
  if (allocated(ibpr)) call mma_deallocate(ibpr)
  if (allocated(ixbpr)) call mma_deallocate(ixbpr)
  if (allocated(iconfs)) call mma_deallocate(iconfs)
  if (allocated(idetvb)) call mma_deallocate(idetvb)
  if (allocated(ia12ind)) call mma_deallocate(ia12ind)
  if (allocated(ib12ind)) call mma_deallocate(ib12ind)
  if (allocated(iapr1)) call mma_deallocate(iapr1)
  if (allocated(ixapr1)) call mma_deallocate(ixapr1)
  if (allocated(ibpr1)) call mma_deallocate(ibpr1)
  if (allocated(ixbpr1)) call mma_deallocate(ixbpr1)
  if (allocated(hessorb)) call mma_deallocate(hessorb)
  if (allocated(hesst)) call mma_deallocate(hesst)
  if (allocated(grad1)) call mma_deallocate(grad1)
  if (allocated(grad2)) call mma_deallocate(grad2)
  if (allocated(gradx)) call mma_deallocate(gradx)
  if (allocated(vec1)) call mma_deallocate(vec1)
  if (allocated(symelm)) call mma_deallocate(symelm)
  if (allocated(iorbrel)) call mma_deallocate(iorbrel)
  if (allocated(north)) call mma_deallocate(north)
  if (allocated(corth)) call mma_deallocate(corth)
  if (allocated(irels)) call mma_deallocate(irels)
  if (allocated(relorb)) call mma_deallocate(relorb)
  if (allocated(ifxorb)) call mma_deallocate(ifxorb)
  if (allocated(ifxstr)) call mma_deallocate(ifxstr)
  if (allocated(idelstr)) call mma_deallocate(idelstr)
  if (allocated(iorts)) call mma_deallocate(iorts)
  if (allocated(irots)) call mma_deallocate(irots)
  if (allocated(izeta)) call mma_deallocate(izeta)
  if (allocated(trprm)) call mma_deallocate(trprm)
  if (allocated(tconstr)) call mma_deallocate(tconstr)
  if (allocated(ipermzeta)) call mma_deallocate(ipermzeta)
  if (allocated(sstruc)) call mma_deallocate(sstruc)
  if (allocated(sstruc2)) call mma_deallocate(sstruc2)
  if (allocated(wdx)) call mma_deallocate(wdx)
  if (allocated(orbs)) call mma_deallocate(orbs)
  if (allocated(cvb)) call mma_deallocate(cvb)
  if (allocated(vbdet)) call mma_deallocate(vbdet)
  if (allocated(aikcof)) call mma_deallocate(aikcof)
  nullify(bikcof)
  if (allocated(cikcof)) call mma_deallocate(cikcof)
  if (allocated(ikcoff)) call mma_deallocate(ikcoff)
  if (allocated(ifnss1)) call mma_deallocate(ifnss1)
  if (allocated(ifnss2)) call mma_deallocate(ifnss2)
  if (allocated(ndetvbs)) call mma_deallocate(ndetvbs)
end subroutine casvb_free

end module casvb_global
