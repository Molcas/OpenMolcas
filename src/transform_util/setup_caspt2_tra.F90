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

subroutine SetUp_CASPT2_Tra(nSym_,nBas_,nOrb_,nIsh_,nAsh_,nFro_,nDel_,CMO_,lthCMO,LuIntM_,LuHlf1_,LuHlf2_,LuHlf3_)

use caspt2_global, only: CMO, LUHLF1, LUHLF2, LUHLF3, LUINTM, NCMO
use caspt2_module, only: nAsh, nBas, nDel, nFro, nIsh, nOrb, nOsh, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym_, nBas_(8), nOrb_(8), nIsh_(8), nAsh_(8), nFro_(8), nDel_(8), lthCMO, LuIntM_
real(kind=wp), target, intent(in) :: CMO_(lthCMO)
integer(kind=iwp), intent(inout) :: LuHlf1_, LuHlf2_, LuHlf3_

!                                                                      *
!***********************************************************************
!                                                                      *
nSym = nSym_
nBas(1:nSym) = nBas_(1:nSym)
nOrb(1:nSym) = nOrb_(1:nSym)
nFro(1:nSym) = nFro_(1:nSym)
nDel(1:nSym) = nDel_(1:nSym)
nAsh(1:nSym) = nAsh_(1:nSym)
nIsh(1:nSym) = nIsh_(1:nSym)
nOsh(1:nSym) = nAsh_(1:nSym)+nIsh_(1:nSym)
!                                                                      *
!***********************************************************************
!                                                                      *
CMO => CMO_
nCMO = lthCMO
!                                                                      *
!***********************************************************************
!                                                                      *
! Open the temporary files here! This is where they get their final
! unit numbers.

call DANAME_MF_wa(LuHlf1_,'LUHLF1')
call DANAME_MF_wa(LuHlf2_,'LUHLF2')
call DANAME_MF_wa(LuHlf3_,'LUHLF3')
LuHlf1 = LuHlf1_
LuHlf2 = LuHlf2_
LuHlf3 = LuHlf3_

! Observe that LuIntM_ should be opened prior to this call!

LuIntM = LuIntM_
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine SetUp_CASPT2_Tra
