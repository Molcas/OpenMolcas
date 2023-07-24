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

use Data_Structures, only: DSBA_Type
use Definitions, only: wp, iwp

implicit none
private

#include "LenIn.fh"

integer(kind=iwp) :: iPL, iPoVec(9), mAdDel(8), mAdFro(8), mAdOcc(8), mAdVir(8), nBas(8), nDel1(8), nDel2(8), nDsto(8), nFro1(8), &
                     nFro2(8), nnB, nTit
real(kind=wp) :: EMP2, Thr_ghs, VECL2
logical(kind=iwp) :: DelGhost, DoCholesky, DoDF, DoLDF
character(len=80) :: Title(10)
type(DSBA_Type) :: Density, DiaA, Mp2Lagr, WDensity
integer(kind=iwp), allocatable :: iDel(:,:), iFro(:,:)
real(kind=wp), allocatable :: CMO(:), EOrb(:)
real(kind=wp), allocatable, target :: EOcc(:), EVir(:)
character(len=LenIn), allocatable :: NamAct(:)
! Unit Numbers & File Names
integer(kind=iwp) :: LuHLF1 = 40, LuHLF2 = 41, LuHLF3 = 42, LuIntA = 43, LuIntM = 44
character(len=8) :: FnHLF1 = 'LUHLF1', FnHLF2 = 'LUHLF2', FnHLF3 = 'LUHLF3', FnIntA = 'ORDINT', FnIntM = 'MOLINT'

public :: CMO, DelGhost, Density, DiaA, DoCholesky, DoDF, DoLDF, EMP2, EOcc, EOrb, EVir, FnHLF1, FnHLF2, FnHLF3, FnIntA, FnIntM, &
          iDel, iFro, iPL, iPoVec, LuHLF1, LuHLF2, LuHLF3, LuIntA, LuIntM, mAdDel, mAdFro, mAdOcc, mAdVir, MBPT2_Clean, MP2Lagr, &
          NamAct, nBas, nDel1, nDel2, nDsto, nFro1, nFro2, nnB, nTit, Title, Thr_ghs, VECL2, WDensity

contains

subroutine MBPT2_Clean()
  use stdalloc, only: mma_deallocate
  if (allocated(CMO)) call mma_deallocate(CMO)
  if (allocated(EOcc)) call mma_deallocate(EOcc)
  if (allocated(EOrb)) call mma_deallocate(EOrb)
  if (allocated(EVir)) call mma_deallocate(EVir)
  if (allocated(iDel)) call mma_deallocate(iDel)
  if (allocated(iFro)) call mma_deallocate(iFro)
  if (allocated(NamAct)) call mma_deallocate(NamAct)
end subroutine MBPT2_Clean

end module MBPT2_Global
