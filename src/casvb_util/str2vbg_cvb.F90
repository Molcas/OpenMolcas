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
!*  STR2VB[CG] : Structure-to-determinant transformation.              *
!*  VB2STR[CG] : Determinant-to-structure transformation.              *
!*                                                                     *
!*  [C] : Transforms coefficients.                                     *
!*  [G] : Transforms a gradient-type quantity.                         *
!*                                                                     *
!*  The difference between [C] and [G] is mainly important when        *
!*  non-orthogonal spin functions are used (Rumer/Project.)            *
!*                                                                     *
!***********************************************************************

subroutine str2vbg_cvb(cvb,cvbdet)

use casvb_global, only: absym, aikcof, i2s_fr, idetvb, kbasiscvb, nalf_fr, nconfion_fr, ndetvb, ndetvb_fr, nel_fr, nfrag, nMs_fr, &
                        nS_fr, nvb, nvb_fr
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: cvb(nvb), cvbdet(ndetvb)
integer(kind=iwp) :: idetvb_add, ifrag, ioffs_cvb, ioffs_cvbdet, kbs

kbs = nint(aikcof(0))
if (kbs /= kbasiscvb) then
  call mkbiks_cvb()
  kbs = kbasiscvb
end if
ioffs_cvb = 1
ioffs_cvbdet = 1
idetvb_add = 1
do ifrag=1,nfrag
  call str2vb2_cvb(aikcof(1:),cvb(ioffs_cvb),cvbdet(ioffs_cvbdet),2,idetvb(idetvb_add),i2s_fr(1,ifrag),nS_fr(ifrag), &
                   nalf_fr(1,ifrag),nMs_fr(ifrag),absym(1),ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,nel_fr(ifrag),nconfion_fr(0,ifrag))
  idetvb_add = idetvb_add+ndetvb_fr(ifrag)
  ioffs_cvb = ioffs_cvb+nvb_fr(ifrag)
  ioffs_cvbdet = ioffs_cvbdet+ndetvb_fr(ifrag)
end do

return

end subroutine str2vbg_cvb
