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

subroutine Cho_Distrib_Vec(Jin,Jfi,iDV,nV)
!
! Unless you know exactly what you are doing,
! do NOT call this routine directly; use Cho_P_Distrib_Vec instead!

use Para_Info, only: MyRank, nProcs
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: Jin, Jfi
integer(kind=iwp), intent(_OUT_) :: iDV(*)
integer(kind=iwp), intent(out) :: nV
integer(kind=iwp) :: i, iNode

nV = 0
do i=Jin,Jfi
  iNode = mod(i-1,nProcs)
  if (myRank == iNode) then
    nV = nV+1
    iDV(nV) = i
  end if
end do

end subroutine Cho_Distrib_Vec
