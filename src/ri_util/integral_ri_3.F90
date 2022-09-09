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

subroutine Integral_RI_3( &
#                        define _CALLING_
#                        include "int_wrout_interface.fh"
                        )
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                    else  integral order within symblk: ijkl

use RICD_Info, only: LDF
use RI_glob, only: iSSOff, nBasSh, klS, nSkal_Valence, nSO, SOShl, ShlSO
use Definitions, only: wp, iwp

implicit none
#include "int_wrout_interface.fh"

#include "macros.fh"
unused_var(MapOrg)
unused_var(iSOSym)
unused_var(nSkal)

!                                                                      *
!***********************************************************************
!                                                                      *
if (LDF) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (mSym == 1) then
    call PLF_LDF_3(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas,lBas,kOp,TInt, &
                   nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO,nSkal_Valence,mSym,iSSOff(0,0,klS))
  else
    call WarningMessage(2,'Not implemented yet!')
    call Abend()
    !call IndSft_RI_3(iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,TInt,nTInt,iTOffs,ShlSO, &
    !                 nBasSh,SOShl,nSO,nSkal_Valence,mSym,iSSOff(:,:,klS))
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (mSym == 1) then
    call PLF_RI_3(AOInt,ijkl,iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,jBas,kBas,lBas,kOp,TInt,nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO, &
                  nSkal_Valence,mSym,iSSOff(0,0,klS))
  else
    call IndSft_RI_3(iCmp,iShell,jBas,kBas,lBas,iAO,iAOst,ijkl,SOInt,nSOint,TInt,nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO, &
                     nSkal_Valence,mSym,iSSOff(:,:,klS))
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Integral_RI_3
