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
!***********************************************************************
!*                                                                     *
!*  VB2CI[CF] := VB to CASSCF transformation.                          *
!*  CI2VB[CG] := CASSCF to VB transformation.                          *
!*  VB2CIAF   := Adds first-order term in CVBDET to CAS CI vector.     *
!*                                                                     *
!*  Transformation between VB vector and CASSCF vector (both defined   *
!*  in determinant basis).                                             *
!*                                                                     *
!*  Normally a simple indexed copy of the coefficients is all that     *
!*  is required. For direct-product wavefunctions, however, a          *
!*  partial contraction is also carried out (relevant for first-       *
!*  order changes and gradient back transformation).                   *
!*                                                                     *
!*  [C] : Transforms coefficients.                                     *
!*  [F] : Transforms a first-order change of the coefficients.         *
!*  [G] : Transforms a gradient-type quantity.                         *
!*                                                                     *
!***********************************************************************

subroutine vb2cic_cvb(cvbdet,civec)

use casvb_global, only: iapr, icnt_ci, iform_ci, ixapr, ndet, ndetvb, nfrag, vbdet
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: cvbdet(ndetvb), civec(0:ndet)
integer(kind=iwp) :: ic, icivec
real(kind=wp) :: dum

icivec = nint(civec(0))
ic = 0

if (iform_ci(icivec) /= 0) then
  write(u6,*) ' Unsupported format in VB2CI :',iform_ci(icivec)
  call abend_cvb()
end if
if (nfrag <= 1) then
  call ci2vb2_cvb(civec(1:),cvbdet,iapr,ixapr,dum,1)
else
  call dpci2vb_cvb(civec(1:),cvbdet,vbdet,ic,dum,1)
end if
icnt_ci(icivec) = 0

return

end subroutine vb2cic_cvb
