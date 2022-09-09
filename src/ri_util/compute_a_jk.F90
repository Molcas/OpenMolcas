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

use ExTerm, only: iMP2prpt, LuAVector
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: jVec, kVec, nVec, iOpt
real(kind=wp) :: Ajk, Fac_ij, Fac_kl
#include "exterm.fh"
integer(kind=iwp) :: iAdrA, lTot
real(kind=wp) :: Ajk_mp2, dum(1)
character(len=*), parameter :: SECNAM = 'Compute_A_jk_mp2'

Ajk = Zero
if (imp2prpt == 2) then
  lTot = 1
  iAdrA = nVec*(kVec-1)+jVec
  call dDaFile(LuAVector(iOpt),2,dum,lTot,iAdrA)
  Ajk_mp2 = dum(1)
  Ajk = Ajk+(Ajk_mp2*Fac_kl*Fac_ij)
end if

return

end subroutine Compute_A_jk_mp2
