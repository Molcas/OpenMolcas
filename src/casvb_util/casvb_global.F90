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

integer(kind=iwp), parameter :: lbuf = 512, mxact_mo = 16, mxdep = 200, mxfield = 500, mxfiles = 40, mxfrag = 10, mxI = 20, &
                                mxirrep_mo = 8, mxMs = 20, mxobj = 100, mxopth = 10, mxprm = 100, mxS = 20, mxstep = 200, &
                                nspinb = 7, nstackrep = 50

integer(kind=iwp) :: i_dep_on_j(mxdep), i2s_fr(mxS,mxfrag), iact_mo(mxact_mo), iaddr, iaddrm(mxfield), iastr_fr(mxfrag), &
                     ibstr_fr(mxfrag), ibuf, ibuffer(lbuf), icase6, icase7, icnt, icode(mxstep), idan(mxfiles), idd(9), ifield, &
                     ifollow, ifsc_fr(mxfrag), igrad, iline = 0, ilv(300), inp, inputmode, ioffs(mxobj+1), iopt2step(0:30), &
                     ioptcode(30), ioptim, ioptstep, ip, ipAnchr, ipdd, ipos(mxstep), ipp, ipp12e, ipp12s, ipp7, iprint, iprm, &
                     iroot, is_set = 0, isaddle, isaddledd, istackrep(nstackrep), iter, iter12e, iter12s, iter7, ivrhs, ix(11), &
                     izbuffer(lbuf), j_dep_on_i(mxdep), joffs(mxobj+1), joptstep, jroot, lenline, loopstep, loopstepmx, &
                     lstprm(mxprm), maxd, mnion_fr(mxfrag), mxdav, mxion_fr(mxfrag), mxit, mxrhs, n_div, nact_mo, &
                     nalf_fr(mxMs,mxfrag), nbas_mo, nbasf_mo(mxirrep_mo), nbasi_mo(mxirrep_mo), nbasisq_mo, &
                     nbassqf_mo(mxirrep_mo), nbassqi_mo(mxirrep_mo), nbet_fr(mxMs,mxfrag), nbuf, ncnt, nconf_fr(mxfrag), &
                     nconfion_fr(0:mxI,mxfrag), nda_fr(mxMs,mxfrag), ndb_fr(mxMs,mxfrag), ndep_ij, ndep_ji, ndetvb_fr(mxfrag), &
                     ndetvb2_fr(mxfrag), nel_fr(mxfrag), nfield, nfieldm, nfold, nfrag, nfrdim, nline = 0, nlold = 0, &
                     nMs_fr(mxfrag), nobj, nopth1(2), nopth2(2), noptim, noptstep, nortiter, nortiterdd, nparm, nroot, &
                     nS_fr(mxfrag), nsym_mo, nvb_fr(mxfrag), nvbr_fr(mxfrag), nvecmx, nvguess, nvrestart, nvrhs, nvtot, nword
real(kind=wp) :: aa1, aa2, alftol, cnrm, cnrmtol, corenrg, cpropt(mxopth), cvbnrm_fr(mxfrag), delopth1(2), delopth2(2), dfx(6), &
                 dfxmin(2), dfxtol, dx(3,6), eigwrngtol, exp12tol, expct, f1, f2, f3, f4, file_id, fxbest, grd(3,6), grdwrngtol, &
                 hh, hhaccfac(5,2), hhkeep, hhmax(2), hhopt(mxopth), hhrejfac(2), hhstart, hhtol(2), oaa2, oaa3, orththr, &
                 orththrdd, ovraa_try, ovrab, ovrab_try, recn, resthr, resthrdd, safety, sgn(6), signtol, singul(3), ww, ww_try, &
                 zzacclim(4,2), zzmax(6), zzmin(6), zzrejmax(2), zzrejmin(2)
character(len=300) :: line
character(len=20) :: form2AD, form2AF, formAD, formAF, formChk1, formChk2, formChk3, formcvp, formE, formMXP1, formMXP2, formMXP3, &
                     formMXP4, formMXP5, formMXP6, formroot, formSymW, formVBWnorm
character(len=8) :: charobj(mxobj)
logical(kind=iwp) :: endwhenclose, follow, have_solved_it, maxize, memdebug, mustdeclare, ndres_ok, release(10), scalesmall(2), &
                     up2date(mxobj)

integer(kind=iwp), parameter :: iunset = -1357924680
character(len=*), parameter :: spinb(nspinb) = ['Kotani      ','Serber      ','Rumer       ','Rumer (LT)  ','Projected   ', &
                                                'Determinants','Determinants'], &
                               spinbkw(nspinb) = ['KOTANI  ','SERBER  ','RUMER   ','LTRUMER ','PROJECT ','DET     ', &
                                                  'DETERM  ']

public :: aa1, aa2, alftol, charobj, cnrm, cnrmtol, corenrg, cpropt, cvbnrm_fr, delopth1, delopth2, dfx, dfxmin, dfxtol, dx, &
          eigwrngtol, endwhenclose, exp12tol, expct, f1, f2, f3, f4, file_id, follow, form2AD, form2AF, formAD, formAF, formChk1, &
          formChk2, formChk3, formcvp, formE, formMXP1, formMXP2, formMXP3, formMXP4, formMXP5, formMXP6, formroot, formSymW, &
          formVBWnorm, fxbest, grd, grdwrngtol, have_solved_it, hh, hhaccfac, hhkeep, hhmax, hhopt, hhrejfac, hhstart, hhtol, &
          i_dep_on_j, i2s_fr, iact_mo, iaddr, iaddrm, iastr_fr, ibstr_fr, ibuf, ibuffer, icase6, icase7, icnt, icode, idan, idd, &
          ifield, ifollow, ifsc_fr, igrad, iline, ilv, inp, inputmode, ioffs, iopt2step, ioptcode, ioptim, ioptstep, ip, ipAnchr, &
          ipdd, ipos, ipp, ipp12e, ipp12s, ipp7, iprint, iprm, iroot, is_set, isaddle, isaddledd, istackrep, iter, iter12e, &
          iter12s, iter7, iunset, ivrhs, ix, izbuffer, j_dep_on_i, joffs, joptstep, jroot, lbuf, lenline, line, loopstep, &
          loopstepmx, lstprm, maxd, maxize, memdebug, mnion_fr, mustdeclare, mxact_mo, mxdav, mxdep, mxfield, mxfiles, mxfrag, &
          mxI, mxion_fr, mxirrep_mo, mxit, mxMs, mxobj, mxopth, mxprm, mxrhs, mxS, mxstep, n_div, nact_mo, nalf_fr, nbas_mo, &
          nbasf_mo, nbasi_mo, nbasisq_mo, nbassqf_mo, nbassqi_mo, nbet_fr, nbuf, ncnt, nconf_fr, nconfion_fr, nda_fr, ndb_fr, &
          ndep_ij, ndep_ji, ndetvb_fr, ndetvb2_fr, ndres_ok, nel_fr, nfield, nfieldm, nfold, nfrag, nfrdim, nline, nlold, nMs_fr, &
          nobj, nopth1, nopth2, noptim, noptstep, nortiter, nortiterdd, nparm, nroot, nS_fr, nspinb, nstackrep, nsym_mo, nvb_fr, &
          nvbr_fr, nvecmx, nvguess, nvrestart, nvrhs, nvtot, nword, oaa2, oaa3, orththr, orththrdd, ovraa_try, ovrab, ovrab_try, &
          recn, release, resthr, resthrdd, safety, scalesmall, sgn, signtol, singul, spinb, spinbkw, up2date, ww, ww_try, &
          zzacclim, zzmax, zzmin, zzrejmax, zzrejmin

end module casvb_global
