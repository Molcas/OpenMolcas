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

subroutine Cho_RI_PutInfo(iPass,iRed)

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
use Cholesky, only: LuRed, LuRed_G, LuRst, LuRst_G
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iPass, iRed
integer(kind=iwp) :: iTmp
logical(kind=iwp) :: doSwap

#ifdef _MOLCAS_MPP_
doSwap = (nProcs > 1) .and. Is_Real_Par()
#else
doSwap = .false.
#endif

if (doSwap) then
  iTmp = LuRed
  LuRed = LuRed_G
  call Cho_PutRed(iPass,iRed) ! save reduced set indices on disk
  LuRed = iTmp
  iTmp = LuRst
  LuRst = LuRst_G
  call Cho_WrRstC(iPass) ! save disk addresses etc. on disk
  LuRst = iTmp
else
  call Cho_PutRed(iPass,iRed) ! save reduced set indices on disk
  call Cho_WrRstC(iPass) ! save disk addresses etc. on disk
end if

return

end subroutine Cho_RI_PutInfo
