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

subroutine cimol2vb_cvb(vec,civec)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: vec(*), civec(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
#include "casinfo_cvb.fh"
integer(kind=iwp) :: icioffs, icivec, istate, istsym_d, isyml, iwr, lcim, nci, ncix(mxirrep)
real(kind=wp) :: fac
integer(kind=iwp), external :: mstackr_cvb

iwr = 0
icivec = nint(civec(1))

if (iwr == 0) call fzero(work(iaddr_ci(icivec)),ndet)

icioffs = 0
do istsym_d=1,nstsym_d
  isyml = istsy_d(istsym_d)
  if (isymv(isyml) /= 1) cycle
  call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),istsy_d(istsym_d))
  nci = ncix(1)
  lcim = mstackr_cvb(nci)
  if (iwr == 0) then
    do istate=1,nstats_d(istsym_d)
      if (abs(weight_d(istate,istsym_d)) > 1.0e-20_wp) then
        call fmove_cvb(vec(1+icioffs),work(lcim),nci)
        icioffs = icioffs+nci
        fac = sqrt(weight_d(istate,istsym_d))
        call mol2vbma_cvb(work(iaddr_ci(icivec)),work(lcim),isyml,fac)
      end if
    end do
  else if (iwr == 1) then
    do istate=1,nstats_d(istsym_d)
      if (abs(weight_d(istate,istsym_d)) > 1.0e-20_wp) then
        call vb2mol_cvb(work(iaddr_ci(icivec)),work(lcim),isyml)
        call fmove_cvb(work(lcim),vec(1+icioffs),nci)
        icioffs = icioffs+nci
      end if
    end do
  end if
  call mfreer_cvb(lcim)
end do

return

end subroutine cimol2vb_cvb
