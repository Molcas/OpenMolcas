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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine R2FIP(CHSPC,NCHSPC,WRK,ipWRK,NUMV,nBasT,iSym,iSkip,irc,JREDC)
! Transform the reduced form to the full form in place

use Cholesky, only: INFVEC, nDimRS, nnBstR
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCHSPC, ipWRK(8), NUMV, nBasT, iSym, iSkip(8)
real(kind=wp), intent(inout) :: CHSPC(NCHSPC), WRK(nBasT*nBasT)
integer(kind=iwp), intent(inout) :: irc, JREDC
integer(kind=iwp) :: ipVecL, iVec, jloc, JREDL, kloc, l_NDIMRS, lscr

l_NDIMRS = size(nDIMRS)
kloc = 0
do iVec=1,NUMV
  if (l_NDIMRS < 1) then
    lscr = NNBSTR(iSym,3)
  else
    JREDL = INFVEC(iVec,2,iSym)
    lscr = nDimRS(iSym,JREDL) !! JRED?
  end if
  kloc = kloc+lscr
end do

ipVecL = 1+kloc !! lscr*(JNUM-1)
jloc = (NUMV-1)*nBasT**2+1
do iVec=NUMV,1,-1
  if (l_NDIMRS < 1) then
    lscr = NNBSTR(iSym,3)
  else
    JREDL = INFVEC(iVec,2,iSym)
    lscr = nDimRS(iSym,JREDL) !! JRED?
  end if
  ipVecL = ipVecL-lscr
  WRK(1:nBasT**2) = Zero
  call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,1,1,1,1,iSym,JREDC,2,ipWRK,WRK,iSkip)
  CHSPC(jloc:jloc+nBasT**2-1) = WRK(1:nBasT**2)
  jloc = jloc-nBasT**2
end do

return

end subroutine R2FIP
