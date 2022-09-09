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

subroutine Integral_RI_2( &
#                        define _CALLING_
#                        include "int_wrout_interface.fh"
                        )
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                   else  integral order within symblk: ijkl

use RI_glob, only: iOffA, SO2Ind
use Definitions, only: wp, iwp

implicit none
#include "int_wrout_interface.fh"

#include "macros.fh"
unused_var(MapOrg)
unused_var(iBas)
unused_var(kBas)
unused_var(IJeqKL)
unused_var(iSOSym)
unused_var(nSkal)
unused_var(itOffs)

if (mSym == 1) then
  call PLF_RI_2(AOInt,ijkl,iCmp(2),iCmp(4),iAO,iAOst,jBas,lBas,kOp,TInt,nTInt,SO2Ind,iOffA,nSOs)
else
  call IndSft_RI_2(iCmp,iShell,jBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,nSOs,TInt,nTInt,SO2Ind,iOffA)
end if

return

end subroutine Integral_RI_2
