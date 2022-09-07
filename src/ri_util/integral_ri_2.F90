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

subroutine Integral_RI_2(iCmp,iShell,MapOrg,iBas,jBas,kBas,lBas,kOp,Shijij,IJeqKL,iAO,iAOst,ijkl,AOInt,SOInt,nSOint,iSOSym,nSkal, &
                         nSOs,TInt,nTInt,itOffs,nSym)
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                   else  integral order within symblk: ijkl

use Wrj12, only: iOffA, SO2Ind
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iCmp(4), iShell(4), MapOrg(4), iBas, jBas, kBas, lBas, kOp(4), iAO(4), iAOst(4), ijkl, nSOint, nSOs, &
                     iSOSym(2,nSOs), nSkal, nTInt, nSym, itOffs(0:nSym-1,0:nSym-1,0:nSym-1)
logical(kind=iwp) :: Shijij, IJeqKL
real(kind=wp) :: AOInt(*), SOInt(*), TInt(nTInt)

if (nSym == 1) then
  call PLF_RI_2(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas,lBas,kOp,TInt, &
                nTInt,SO2Ind,iOffA,nSOs)
else
  call IndSft_RI_2(iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,TInt,nTInt,iTOffs,SO2Ind,iOffA)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(MapOrg)
  call Unused_integer(nSkal)
end if

end subroutine Integral_RI_2
