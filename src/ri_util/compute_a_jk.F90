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
! Copyright (C) Jonas Bostrom                                          *
!***********************************************************************

subroutine Compute_A_jk_mp2(jVec,kVec,Ajk,fac_ij,fac_kl,nVec,iOpt)
!***********************************************************************
!   Author: J Bostrom                                                  *
!                                                                      *
!   Purpose: Loading A-matrix for mp2 from disk                        *
!                                                                      *
!***********************************************************************

use RI_glob, only: iMP2prpt, LuAVector
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: jVec, kVec, nVec, iOpt
real(kind=wp), intent(out) :: Ajk
real(kind=wp), intent(in) :: Fac_ij, Fac_kl
integer(kind=iwp) :: iAdrA
real(kind=wp) :: Ajk_mp2(1)
character(len=*), parameter :: SECNAM = 'Compute_A_jk_mp2'

if (iMP2prpt == 2) then
  iAdrA = nVec*(kVec-1)+jVec
  call dDaFile(LuAVector(iOpt),2,Ajk_mp2,1,iAdrA)
  Ajk = Ajk_mp2(1)*Fac_kl*Fac_ij
else
  Ajk = Zero
end if

return

end subroutine Compute_A_jk_mp2
