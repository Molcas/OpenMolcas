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

module MBPT2_Global

use Definitions, only: wp, iwp

implicit none
private

#include "LenIn.fh"

integer(kind=iwp) :: ip_Density(8), ip_DiaA(8), ip_First_Density, ip_First_DiaA, ip_First_Mp2Lagr, ip_First_WDensity, &
                     ip_Mp2Lagr(8), ip_WDensity(8), ipCMO, ipEOcc, ipEVir, ipInt1, ipInt1_2, ipInt2, ipInt2_2, iPL, iPoVec(9), &
                     ipScr1, l_Density, l_DiaA, l_Mp2Lagr, mAdDel(8), mAdFro(8), mAdOcc(8), mAdVir(8), nBas(8), nDel1(8), &
                     nDel2(8), nDsto(8), nFro1(8), nFro2(8), nnB, nTit
real(kind=wp) :: EMP2, Thr_ghs, VECL2
logical(kind=iwp) :: DelGhost, DoCholesky, DoDF, DoLDF
character(len=80) :: Title(10)
integer(kind=iwp), allocatable :: iDel(:,:), iFro(:,:)
character(len=LenIn), allocatable :: NamAct(:)
! Unit Numbers & File Names
integer(kind=iwp) :: LuHLF1 = 40, LuHLF2 = 41, LuHLF3 = 42, LuIntA = 43, LuIntM = 44
character(len=8) :: FnHLF1 = 'LUHLF1', FnHLF2 = 'LUHLF2', FnHLF3 = 'LUHLF3', FnIntA = 'ORDINT', FnIntM = 'MOLINT'

public :: EMP2, DelGhost, DoCholesky, DoDF, DoLDF, FnHLF1, FnHLF2, FnHLF3, FnIntA, FnIntM, iDel, iFro, ip_Density, ip_DiaA, &
          ip_First_Density, ip_First_DiaA, ip_First_Mp2Lagr, ip_First_WDensity, ip_Mp2Lagr, ip_WDensity, ipCMO, ipEOcc, ipEVir, &
          ipInt1, ipInt1_2, ipInt2, ipInt2_2, iPL, iPoVec, ipScr1, l_Density, l_DiaA, l_Mp2Lagr, LuHLF1, LuHLF2, LuHLF3, LuIntA, &
          LuIntM, mAdDel, mAdFro, mAdOcc, mAdVir, MBPT2_Clean, NamAct, nBas, nDel1, nDel2, nDsto, nFro1, nFro2, nnB, nTit, Title, &
          Thr_ghs, VECL2

contains

subroutine MBPT2_Clean()
  use stdalloc, only: mma_deallocate
  if (allocated(iDel)) call mma_deallocate(iDel)
  if (allocated(iFro)) call mma_deallocate(iFro)
  if (allocated(NamAct)) call mma_deallocate(NamAct)
end subroutine MBPT2_Clean

end module MBPT2_Global
