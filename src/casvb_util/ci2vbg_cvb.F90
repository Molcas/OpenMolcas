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

subroutine ci2vbg_cvb(civec,cvbdet)

use casvb_global, only: nfrag

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension cvbdet(ndetvb), civec(*)

icivec = nint(civec(1))
ic = 2

if (iform_ci(icivec) /= 0) then
  write(6,*) ' Unsupported format in CI2VB :',iform_ci(icivec)
  call abend_cvb()
end if
if (nfrag <= 1) then
  call ci2vb2_cvb(work(iaddr_ci(icivec)),cvbdet,iwork(ll(11)),iwork(ll(12)),dum,0)
else
  call dpci2vb_cvb(work(iaddr_ci(icivec)),cvbdet,work(lv(5)),ic,dum,0)
end if

return

end subroutine ci2vbg_cvb
