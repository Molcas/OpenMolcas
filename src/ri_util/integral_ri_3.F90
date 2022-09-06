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

subroutine Integral_RI_3(iCmp,iShell,MapOrg,iBas,jBas,kBas,lBas,kOp,Shijij,IJeqKL,iAO,iAOst,ijkl,AOInt,SOInt,nSOint,iSOSym,nSkal, &
                         nSOs,TInt,nTInt,itOffs,nSym)
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                    else  integral order within symblk: ijkl

use RICD_Info, only: LDF
use j12

implicit real*8(A-H,O-Z)
real*8 AOInt(*), SOInt(*), TInt(nTInt)
integer iCmp(4), iShell(4), iAO(4), iAOst(4), kOp(4), iSOSym(2,nSOs), itOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4)
logical Shijij, IJeqKL

!                                                                      *
!***********************************************************************
!                                                                      *
if (LDF) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (nSym == 1) then
    call PLF_LDF_3(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas,lBas,kOp,TInt, &
                   nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO,nSkal_Valence,nSym,iSSOff(0,0,klS))
  else
    call WarningMessage(2,'Not implemented yet!')
    call Abend()
    !call IndSft_RI_3(iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,TInt,nTInt,iTOffs,ShlSO, &
    !                 nBasSh,SOShl,nSO,nSkal_Valence,nSym,iSSOff(:,:,klS))
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (nSym == 1) then
    call PLF_RI_3(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas,lBas,kOp,TInt, &
                  nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO,nSkal_Valence,nSym,iSSOff(0,0,klS))
  else
    call IndSft_RI_3(iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,TInt,nTInt,iTOffs,ShlSO, &
                     nBasSh,SOShl,nSO,nSkal_Valence,nSym,iSSOff(:,:,klS))
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(MapOrg)
  call Unused_integer(nSkal)
end if

end subroutine Integral_RI_3
