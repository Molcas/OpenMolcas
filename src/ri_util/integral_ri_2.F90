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
!                   else   integral order within symblk: ijkl

use RI_glob, only: iOffA, SO2Ind
use Definitions, only: wp, iwp

implicit none
#include "int_wrout_interface.fh"
integer(kind=iwp) :: iAO(4), iAOst(4), iBas, iCmp(4), iShell(4), jBas, kBas, kOp(4), lBas
logical(kind=iwp) :: Shijij

#include "macros.fh"
unused_var(iBas)
unused_var(kBas)
unused_var(iSOSym)

iCmp(:) = iSD4(2,:)
iShell(:) = iSD4(11,:)
iAO(:) = iSD4(7,:)
iAOst(:) = iSD4(8,:)
iBas = iSD4(19,1)
jBas = iSD4(19,2)
kBas = iSD4(19,3)
lBas = iSD4(19,4)
Shijij = (iSD4(0,1) == iSD4(0,3)) .and. (iSD4(10,1) == iSD4(10,3)) .and. (iSD4(0,2) == iSD4(0,4)) .and. (iSD4(10,2) == iSD4(10,4))

if (mSym == 1) then
  kOp(:) = 0
  call PLF_RI_2(AOInt,ijkl,iCmp(2),iCmp(4),iAO,iAOst,jBas,lBas,kOp,TInt,nTInt,SO2Ind,iOffA,nSOs)
else
  call IndSft_RI_2(iCmp,iShell,jBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint,nSOs,TInt,nTInt,SO2Ind,iOffA)
end if

end subroutine Integral_RI_2
