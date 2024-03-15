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

! What calculation? :
!--------------------
! endvar nmcscf service variat

! Definition of CASSCF active space :
!------------------------------------
! ... Orbital space ...
! ityp
! ... Orbital-space symmetries ...
! isym isympr isymv nirrep nsym
! ... State-averaged wavefunction ...
! esym
! ... Number of determinants in each irrep ...
! ncivb

! Numbers defining, or relating to, the active CI space :
!--------------------------------------------------------
! n1a n1b nalf nam1 nbet nbm1 nda ndb ndet nel noe norb

! ... Utilize alpha<->beta symmetry in various circumstances ...
! absym

! Numbers relating to the definition of the VB wavefunction :
!------------------------------------------------------------
! kbasis kbasiscvb mnion mxion naprodvb nbprodvb nconf ndetvb npvb nvb nvbinp sc

! CASVB algorithm (for evaluating g/G, etc.) :
!---------------------------------------------
! ... Overlap-based or energy based? ...
! icrit
! ... Projection operators ...
! proj projcas projsym

! Parameters relating to the optimization procedure :
!----------------------------------------------------
! imethod isaddle maxdav mxiter

! Optimization handler :
!-----------------------
! convinone ifinish initial ioptc_new

! Orbital permutation :
!----------------------
! iorbprm

! CI vectors :
!-------------
! icnt_ci iform_ci mxciobj ndres nv

! Symmetry information and constraints:
!--------------------------------------
! ... General ...
! sym
! ... Symmetry elements ...
! mxsyme nsyme tags
! ... Symmetry-constrained orbitals ...
! mxops ndimrel nijrel niorth norbrel
! ... Fixed orbitals ...
! nfxorb
! ... Orthogonalized orbitals/deleted rotations ...
! ndrot nort
! ... Symmetry-constrained/fixed/deleted structures ...
! iconstruc lfxvb lzrvb nconstr nfxvb nzeta nzrvb

! Various maximum dimensions defined in input :
!----------------------------------------------
! mxnvb

! All2free/free2all variable transformation :
!--------------------------------------------
! nfr nfrorb nfrvb npr nprorb nprvb orbfr_is_unit orbopt strucopt

! Integrals (number of, core energy, etc.) :
!-------------------------------------------
! mxaobf

! Analysis :
!-----------
! iciweights ishstruc ivbweights lcalccivbs lcalcevb lcalcsvb lciweights npcf sij

! Strictly localized calculations :
!----------------------------------
! plc_const ploc

! Various quantities calculated from VB wavefunction :
!-----------------------------------------------------
! cvbnrm evb ovraa svb

! Usage statistics :
!-------------------
! cpu0, cpu_prev n_2el n_applyh n_applyt n_cihess n_hess n_iter n_orbhess

! Print control :
!----------------
! ipr iprec iwidth

! Save memory? :
!---------------
! dxmove memplenty

! Format statements depending on iprec:
!--------------------------------------
! form2AD form2AF formAD formAF formChk1 formChk2 formChk3 formcvp formE formMXP1 formMXP2 formMXP3 formMXP4 formMXP5 formMXP6
! formroot formSymW formVBWnorm

! _c quantities are based on latest CASSCF
! _d quantities are those actually used (= _c if not changed in input)
! iorclos_* iorcore_* iorocc_* istms2_* istnel_* istsy_* mcore_* nstats_* nstsym_* weight_*

use Definitions, only: wp, iwp

implicit none
private

type gjorb_type
  real(kind=wp), allocatable :: r(:,:)
  integer(kind=iwp), allocatable :: i1(:)
  integer(kind=iwp), allocatable :: i2(:,:)
end type gjorb_type

integer(kind=iwp), parameter :: lbuf = 512, max_rec = 5000, mxact_mo = 16, mxciobj = 20, mxdep = 200, mxfield = 500, mxfiles = 40, &
                                mxfrag = 10, mxI = 20, mxirrep = 8, mxirrep_ci = 8, mxirrep_mo = 8, mxMs = 20, mxobj = 100, &
                                mxopth = 10, mxops = 32, mxorb_cvb = 50, mxprm = 100, mxS = 20, mxstep = 200, mxstsy_ci = 8, &
                                mxstt_ci = 20, mxsyme = 32, mxunits = 8, nspinb = 7, nstackrep = 50

integer(kind=iwp) :: i_dep_on_j(mxdep), i2s_fr(mxS,mxfrag), iact_mo(mxact_mo), iaddrm(mxfield), ibuf, ibuffer(lbuf), icase6, &
                     icase7, iciweights, icnt, icnt_ci(mxciobj), icode(mxstep), iconstruc, icrit, idan(mxfiles), ifield, &
                     ifilio(max_rec), ifinish, ifollow, iform_ci(mxciobj), ifsc_fr(mxfrag), ifvb, iline = 0, ilv(300), imethod, &
                     initial, inp, inputmode, invec_cvb, ioffs(mxobj+1), iopt2step(0:30), ioptc_new, ioptcode(30), ioptim, &
                     ioptstep, iorbprm(mxorb_cvb), iorclos_c(mxirrep_ci), iorclos_d(mxirrep_ci), iorcore_c(mxirrep_ci), &
                     iorcore_d(mxirrep_ci), iorder(mxunits), iorocc_c(mxirrep_ci), iorocc_d(mxirrep_ci), ip, ipAnchr, ipdd, &
                     ipos(mxstep), ipp10, ipp12e, ipp12s, ipp7, ipr(10), iprec, iprint, iprm, iroot, is_set = 0, &
                     isaddle, isaddledd, isaddleo, ishstruc, istackrep(nstackrep), istms2_c(mxstsy_ci), istms2_d(mxstsy_ci), &
                     istnel_c(mxstsy_ci), istnel_d(mxstsy_ci), istsy_c(mxstsy_ci), istsy_d(mxstsy_ci), isym, isympr(mxirrep), &
                     isymv(mxirrep), iter10, iter12e, iter12s, iter7, ityp(mxorb_cvb), ivbweights, iwidth, izbuffer(lbuf), &
                     j_dep_on_i(mxdep), joffs(mxobj+1), joptstep, jroot, kbasis, kbasiscvb, &
                     lenline, lfxvb, loopstep, &
                     loopstepmx, lstprm(mxprm), lzrvb, maxd, maxdav, mcore_c, mcore_d, &
                     mnion, mnion_fr(mxfrag), mxaobf, mxdav, mxion, mxion_fr(mxfrag), mxit, mxiter, mxnvb = 0, mxrhs, n1a, n1b, &
                     n_2el, n_applyh, n_applyt, n_cihess, n_div, n_hess, n_iter, n_orbhess, nact_mo, nalf, nalf_fr(mxMs,mxfrag), &
                     nam1, naprodvb, nbas_mo, nbasf_mo(mxirrep_mo), nbasi_mo(mxirrep_mo), nbasisq_mo, nbassqf_mo(mxirrep_mo), &
                     nbassqi_mo(mxirrep_mo), nbet, nbet_fr(mxMs,mxfrag), nbm1, nbprodvb, nbuf, ncivb(mxirrep), ncnt, nconf, &
                     nconf_fr(mxfrag), nconfion_fr(0:mxI,mxfrag), nconstr, nda, nda_fr(mxMs,mxfrag), ndb, ndb_fr(mxMs,mxfrag), &
                     ndep_ij, ndep_ji, ndet, ndetvb, ndetvb_fr(mxfrag), ndetvb2_fr(mxfrag), ndimrel, ndres, ndrot, nel, &
                     nel_fr(mxfrag), nfield, nfold, nfr, nfrag, nfrdim, nfrorb, nfrvb, nfxorb, nfxvb, nijrel, niorth, nirrep, &
                     nline = 0, nlold = 0, nmcscf, nMs_fr(mxfrag), nobj, noe, nopth1(2), nopth2(2), noptim, noptstep, norb, &
                     norbrel, nort, nortiter, nortiterdd, nparm, npcf, npr, nprorb, nprvb, npvb, nrec, nroot, nS_fr(mxfrag), &
                     nstats_c(mxstsy_ci), nstats_d(mxstsy_ci), nstsym_c, nstsym_d, nsym, nsym_mo, nsyme, nv, nvb, nvb_fr(mxfrag), &
                     nvbinp, nvbr_fr(mxfrag), nvecmx, nvguess, nvrestart, nvrhs, nvtot, nword, nzeta, nzrvb
real(kind=wp) :: aa1, aa2, alftol, cnrm, cnrmtol, corenrg, cpropt(mxopth), cpu0, cpu_prev, cvbnrm, cvbnrm_fr(mxfrag), delopth1(2), &
                 delopth2(2), dfx(6), dfxmin(2), dfxtol, dx(3,6), eigwrngtol, esym(mxirrep), evb, exp12tol, expct, f1, f2, f3, f4, &
                 file_id, fileids(max_rec), fxbest, grd(3,6), grdwrngtol, hh, hhaccfac(5,2), hhkeep, hhmax(2), hhopt(mxopth), &
                 hhrejfac(2), hhstart, hhtol(2), oaa2, oaa3, orththr, orththrdd, ovraa, ovraa_try, ovrab, ovrab_try, recinp, &
                 recinp_old, recn, recn_jobiph, recn_jobold, recn_oneint, recn_tmp01, recn_tmp02, recn_tmp03, recn_tmp04, &
                 recn_vbwfn, resthr, resthrdd, safety, savvb, savvbci, sgn(6), signtol, singul(3), strtci, strtint, strtmo, &
                 strtvb, svb, thresh_io, weight_c(mxstt_ci,mxstsy_ci), weight_d(mxstt_ci,mxstsy_ci), ww, ww_try, zzacclim(4,2), &
                 zzmax(6), zzmin(6), zzrejmax(2), zzrejmin(2)
logical(kind=iwp) :: absym(5), convinone, dxmove, endvar, endwhenclose, follow, have_solved_it, lcalccivbs, lcalcevb, lcalcsvb, &
                     lciweights, maxize, memplenty, mustdeclare, ndres_ok, orbfr_is_unit, orbopt, plc_const, ploc, proj, projcas, &
                     projsym, release(10), sc, scalesmall(2), service, sij, strucopt, sym, up2date(mxobj), variat
character(len=300) :: line
character(len=20) :: filename(max_rec), form2AD, form2AF, formAD, formAF, formChk1, formChk2, formChk3, formcvp, formE, formMXP1, &
                     formMXP2, formMXP3, formMXP4, formMXP5, formMXP6, formroot, formSymW, formVBWnorm
character(len=8) :: charobj(mxobj)
character(len=3) :: tags(mxsyme)
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
                          civb5(:) => null(), civb6(:) => null(), civb7(:) => null(), civb8(:) => null(), phbto(:,:) => null()

integer(kind=iwp), parameter :: iunset = -1357924680
logical(kind=iwp), parameter :: ifhamil = .true., ifmos = .true.
character(len=*), parameter :: spinb(nspinb) = ['Kotani      ','Serber      ','Rumer       ','Rumer (LT)  ','Projected   ', &
                                                'Determinants','Determinants'], &
                               spinbkw(nspinb) = ['KOTANI  ','SERBER  ','RUMER   ','LTRUMER ','PROJECT ','DET     ', &
                                                  'DETERM  ']

public :: aa1, aa2, absym, aikcof, alftol, ap, axc, bikcof, c, casvb_free, charobj, cikcof, civb1, civb2, civb3, civb4, civb5, &
          civb6, civb7, civb8, civbvecs, cnrm, cnrmtol, confsinp, convinone, corenrg, corth, cpropt, cpu0, cpu_prev, cvb, cvbdet, &
          cvbnrm, cvbnrm_fr, cvbsspn, cvbstot, cvbtry, delopth1, delopth2, dfx, dfxmin, dfxtol, dvbdet, dx, dxmove, eigval, &
          eigvec, eigwrngtol, endvar, endwhenclose, esym, evb, evbdet, exp12tol, expct, f1, f2, f3, f4, file_id, fileids, &
          filename, follow, form2AD, form2AF, formAD, formAF, formChk1, formChk2, formChk3, formcvp, formE, formMXP1, formMXP2, &
          formMXP3, formMXP4, formMXP5, formMXP6, formroot, formSymW, formVBWnorm, fxbest, gjorb, gjorb_type, gjorb2, gjorb3, &
          grad1, grad2, gradx, grd, grdwrngtol, gsinp, have_solved_it, hessorb, hesst, hh, hhaccfac, hhkeep, hhmax, hhopt, &
          hhrejfac, hhstart, hhtol, i_dep_on_j, i1alf, i1bet, i1c, i2s_fr, ia12ind, iact_mo, iaddrm, iafrm, iapr, iapr1, iato, &
          ib12ind, ibfrm, ibpr, ibpr1, ibto, ibuf, ibuffer, icase6, icase7, icfrm, iciweights, icnt, icnt_ci, icode, iconfs, &
          iconstruc, icrit, icto, idan, idelstr, idetvb, ifhamil, ifield, ifilio, ifinish, ifmos, ifnss1, ifnss2, ifollow, &
          iform_ci, ifsc_fr, ifvb, ifxorb, ifxstr, ikcoff, iline, ilv, imethod, initial, inp, inputmode, invec_cvb, ioffs, &
          iopt2step, ioptc_new, ioptcode, ioptim, ioptstep, iorbprm, iorbrel, iorclos_c, iorclos_d, iorcore_c, iorcore_d, iorder, &
          iorocc_c, iorocc_d, iorts, ip, ipAnchr, ipdd, ipermzeta, ipos, ipp10, ipp12e, ipp12s, ipp7, ipr, iprec, &
          iprint, iprm, irels, iroot, irots, is_set, isaddle, isaddledd, isaddleo, ishstruc, istackrep, istms2_c, istms2_d, &
          istnel_c, istnel_d, istsy_c, istsy_d, isym, isympr, isymv, iter10, iter12e, iter12s, iter7, ityp, iunset, ivbweights, &
          iwidth, ixapr, ixapr1, ixbpr, ixbpr1, izbuffer, izeta, j_dep_on_i, joffs, joptstep, jroot, kbasis, kbasiscvb, lbuf, &
          lcalccivbs, lcalcevb, lcalcsvb, lciweights, &
          lenline, lfxvb, line, loopstep, loopstepmx, lstprm, &
          lzrvb, max_rec, maxd, maxdav, maxize, mcore_c, mcore_d, memplenty, mnion, mnion_fr, mustdeclare, mxact_mo, mxaobf, &
          mxdav, mxdep, mxfield, mxfiles, mxfrag, mxI, mxion, mxion_fr, mxirrep, mxirrep_ci, mxirrep_mo, mxit, mxiter, mxMs, &
          mxnvb, mxobj, mxopth, mxops, mxorb_cvb, mxprm, mxrhs, mxS, mxstep, mxstsy_ci, mxstt_ci, mxsyme, mxunits, n1a, n1b, &
          n_2el, n_applyh, n_applyt, n_cihess, n_div, n_hess, n_iter, n_orbhess, nact_mo, nalf, nalf_fr, nam1, naprodvb, nbas_mo, &
          nbasf_mo, nbasi_mo, nbasisq_mo, nbassqf_mo, nbassqi_mo, nbet, nbet_fr, nbm1, nbprodvb, nbuf, ncivb, ncnt, nconf, &
          nconf_fr, nconfion_fr, nconstr, nda, nda_fr, ndb, ndb_fr, ndep_ij, ndep_ji, ndet, ndetvb, ndetvb_fr, ndetvb2_fr, &
          ndetvbs, ndimrel, ndres, ndres_ok, ndrot, nel, nel_fr, nfield, nfold, nfr, nfrag, nfrdim, nfrorb, nfrvb, nfxorb, nfxvb, &
          nijrel, niorth, nirrep, nline, nlold, nmcscf, nMs_fr, nobj, noe, nopth1, nopth2, noptim, noptstep, norb, norbrel, nort, &
          north, nortiter, nortiterdd, nparm, npcf, npr, nprorb, nprvb, npvb, nrec, nroot, nS_fr, nspinb, nstackrep, nstats_c, &
          nstats_d, nstsym_c, nstsym_d, nsym, nsym_mo, nsyme, nv, nvb, nvb_fr, nvbinp, nvbr_fr, nvecmx, nvguess, nvrestart, nvrhs, &
          nvtot, nword, nzeta, nzrvb, oaa2, oaa3, odx, odxp, ograd, ogradp, orbfr_is_unit, orbinv, orbopt, orbs, orbstry, orththr, &
          orththrdd, ovraa, ovraa_try, ovrab, ovrab_try, owrk, owrk2, phato, phbto, phcto, plc_const, ploc, proj, projcas, &
          projsym, recinp, recinp_old, recn, recn_jobiph, recn_jobold, recn_oneint, recn_tmp01, recn_tmp02, recn_tmp03, &
          recn_tmp04, recn_vbwfn, release, relorb, res, resthr, resthrdd, rhs, rhsp, savvb, savvbci, safety, sc, scalesmall, &
          service, sgn, signtol, sij, singul, solp, solp_res, sorbs, span, spinb, spinbkw, sstruc, sstruc2, strtci, strtint, &
          strtmo, strtvb, strucopt, svb, sxc, sym, symelm, tags, tconstr, thresh_io, trprm, up2date, variat, vbdet, vec1, wdx, &
          weight_c, weight_d, ww, ww_try, zzacclim, zzmax, zzmin, zzrejmax, zzrejmin

contains

subroutine casvb_free()
  use stdalloc, only: mma_deallocate
  if (allocated(civbvecs)) call mma_deallocate(civbvecs)
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
