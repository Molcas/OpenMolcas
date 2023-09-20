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
!**************************************************************
!** Transformation between VB structures and VB determinants **
!** Also generation of spin-function to determinant transf.  **
!** and of perfect-pairing guess.                            **
!**************************************************************
!***********************************************************************
!*                                                                     *
!*  Transformation between VB structures (VB CSFs) and determinants.   *
!*                                                                     *
!*  The VB structures are ordered according to spatial                 *
!*  configurations, and are normally defined in terms of spin          *
!*  eigenfunctions.                                                    *
!*  VB determinants are ordered with alpha & beta indices in           *
!*  increasing order, with alpha being the slower index.               *
!*                                                                     *
!*  STR2VB[CFG] : Structure-to-determinant transformation.             *
!*  VB2STR[CG] : Determinant-to-structure transformation.              *
!*                                                                     *
!*  [C] : Transforms coefficients.                                     *
!*  [F] : Transforms a first-order change of the coefficients.         *
!*  [G] : Transforms a gradient-type quantity.                         *
!*                                                                     *
!*  At present there is no actual difference between [C] and [F].      *
!*                                                                     *
!*  The difference between [CF] and [G] is mainly important when       *
!*  non-orthogonal spin functions are used (Rumer/Project.)            *
!*                                                                     *
!***********************************************************************

subroutine str2vbg_cvb(cvb,cvbdet)

use casvb_global, only: i2s_fr, mnion_fr, mxion_fr, nalf_fr, nconf_fr, nconfion_fr, ndetvb_fr, nel_fr, nfrag, nMs_fr, nS_fr, nvb_fr

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension cvb(nvb), cvbdet(ndetvb)

kab = 1

kbs = nint(work(lb(kab)))
if (kbs /= kbasiscvb) then
  call mkbiks_cvb()
  kbs = kbasiscvb
end if
ioffs_cvb = 1
ioffs_cvbdet = 1
idetvb_add = ll(17)
ifnss_add = lb(4)
if (kbasiscvb == 6) ifnss_add = lb(5)
ndetvbs_add = lb(6)
do ifrag=1,nfrag
  iwrk = mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
  call str2vb2_cvb(work(lb(kab)+1),iwork(lb(3)),cvb(ioffs_cvb),cvbdet(ioffs_cvbdet),2,iwork(idetvb_add),i2s_fr(1,ifrag), &
                   nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),iwork(ifnss_add),iwork(ndetvbs_add),absym(1),mnion_fr(ifrag), &
                   mxion_fr(ifrag),nconf_fr(ifrag),ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,nel_fr(ifrag),nalf_fr(1,ifrag),nel, &
                   work(iwrk),nconfion_fr(0,ifrag))
  call mfreer_cvb(iwrk)
  idetvb_add = idetvb_add+ndetvb_fr(ifrag)
  ioffs_cvb = ioffs_cvb+nvb_fr(ifrag)
  ioffs_cvbdet = ioffs_cvbdet+ndetvb_fr(ifrag)
end do

return

end subroutine str2vbg_cvb
