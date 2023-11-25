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

! This subroutine should be in a module
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Integral_RI_3( &
#                        define _CALLING_
#                        include "int_wrout_interface.fh"
                        )
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                    else  integral order within symblk: ijkl

use RI_glob, only: iSSOff, nBasSh, klS, nSkal_Valence, nSO, SOShl, ShlSO
use Int_Options, only: iTOffs

implicit none
#include "int_wrout_interface.fh"

#include "macros.fh"
unused_var(Shijij)
unused_var(iSOSym)
unused_var(iBas)

!                                                                      *
!***********************************************************************
!                                                                      *
if (mSym == 1) then
  call PLF_RI_3(AOInt,ijkl,iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,jBas,kBas,lBas,kOp,TInt,nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO, &
                nSkal_Valence,mSym,iSSOff(0,0,klS))
else
  call IndSft_RI_3(iCmp,iShell,jBas,kBas,lBas,iAO,iAOst,ijkl,SOInt,nSOint,TInt,nTInt,iTOffs,ShlSO,nBasSh,SOShl,nSO,nSkal_Valence, &
                   mSym,iSSOff(:,:,klS))
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Integral_RI_3
