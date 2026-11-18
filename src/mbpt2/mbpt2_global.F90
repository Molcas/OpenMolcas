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
use Molcas, only: LenIn
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iPL, iPoVec(9), mAdDel(8), mAdFro(8), mAdOcc(8), mAdVir(8), nBas(8), nDel1(8), nDel2(8), nDsto(8), nFro1(8), &
                     nFro2(8), nnB, nTit
real(kind=wp) :: EMP2, Thr_ghs, VECL2
logical(kind=iwp) :: DelGhost, DoCholesky, DoDF
character(len=80) :: Title(10)
type(DSBA_Type) :: Density, DiaA, Mp2Lagr, WDensity
integer(kind=iwp), allocatable :: iDel(:,:), iFro(:,:)
real(kind=wp), allocatable :: EOrb(:)
real(kind=wp), allocatable, target :: CMO_Internal(:), EOcc(:), EVir(:)
real(kind=wp), pointer:: CMO(:) => null()
character(len=LenIn), allocatable :: NamAct(:)
! Unit Numbers & File Names
integer(kind=iwp) :: LuHLF1 = 40, LuHLF2 = 41, LuHLF3 = 42, LuIntA = 43, LuIntM = 44
character(len=8) :: FnHLF1 = 'LUHLF1', FnHLF2 = 'LUHLF2', FnHLF3 = 'LUHLF3', FnIntA = 'ORDINT', FnIntM = 'MOLINT'

public :: CMO, CMO_Internal, DelGhost, Density, DiaA, DoCholesky, DoDF, EMP2, EOcc, EOrb, EVir, FnHLF1, FnHLF2, FnHLF3, FnIntA, &
          FnIntM, iDel, iFro, iPL, iPoVec, LuHLF1, LuHLF2, LuHLF3, LuIntA, LuIntM, mAdDel, mAdFro, mAdOcc, mAdVir, MBPT2_Clean, &
          MP2Lagr, MBPT2_All_Finite, NamAct, nBas, nDel1, nDel2, nDsto, nFro1, nFro2, nnB, nTit, Title, Thr_ghs, VECL2, WDensity

interface MBPT2_All_Finite
  module procedure MBPT2_All_Finite_Scalar
  module procedure MBPT2_All_Finite_Vector
  module procedure MBPT2_All_Finite_Matrix
  module procedure MBPT2_All_Finite_Tensor
end interface MBPT2_All_Finite

contains


! Shared finite-value check for UMP2/ROMP2.  Chk4NaN follows the
! established OpenMolcas call pattern (n, one-dimensional data, ierr).
! Explicit range comparisons additionally reject positive/negative infinity,
! even if the installed Chk4NaN implementation detects only NaNs.
! Unlike the IEEE intrinsic module, this works with pre-GCC-5 builds.
logical function MBPT2_All_Finite_Scalar(X) result(IsFinite)
  real(kind=wp), intent(in) :: X
  real(kind=wp) :: Values(1)
  integer(kind=iwp) :: iErr
  Values(1)=X
  call Chk4NaN(1,Values,iErr)
  IsFinite=(iErr==0)
  if (IsFinite) IsFinite=(X<=huge(X)).and.(X>=-huge(X))
end function MBPT2_All_Finite_Scalar

logical function MBPT2_All_Finite_Vector(X) result(IsFinite)
  real(kind=wp), intent(in) :: X(:)
  integer(kind=iwp) :: iErr
  real(kind=wp) :: Values(size(X))
  Values(:)=X(:)
  iErr=0
  if (size(Values)>0) call Chk4NaN(size(Values),Values,iErr)
  IsFinite=(iErr==0)
  if (IsFinite) IsFinite=all((X<=huge(0.0_wp)).and.(X>=-huge(0.0_wp)))
end function MBPT2_All_Finite_Vector

logical function MBPT2_All_Finite_Matrix(X) result(IsFinite)
  real(kind=wp), intent(in) :: X(:,:)
  integer(kind=iwp) :: iErr
  real(kind=wp) :: Values(size(X))
  Values(:)=reshape(X,[size(X)])
  iErr=0
  if (size(Values)>0) call Chk4NaN(size(Values),Values,iErr)
  IsFinite=(iErr==0)
  if (IsFinite) IsFinite=all((X<=huge(0.0_wp)).and.(X>=-huge(0.0_wp)))
end function MBPT2_All_Finite_Matrix

logical function MBPT2_All_Finite_Tensor(X) result(IsFinite)
  real(kind=wp), intent(in) :: X(:,:,:)
  integer(kind=iwp) :: iErr
  real(kind=wp) :: Values(size(X))
  Values(:)=reshape(X,[size(X)])
  iErr=0
  if (size(Values)>0) call Chk4NaN(size(Values),Values,iErr)
  IsFinite=(iErr==0)
  if (IsFinite) IsFinite=all((X<=huge(0.0_wp)).and.(X>=-huge(0.0_wp)))
end function MBPT2_All_Finite_Tensor

subroutine MBPT2_Clean()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(CMO_Internal,safe='*')
  call mma_deallocate(EOcc,safe='*')
  call mma_deallocate(EOrb,safe='*')
  call mma_deallocate(EVir,safe='*')
  call mma_deallocate(iDel,safe='*')
  call mma_deallocate(iFro,safe='*')
  call mma_deallocate(NamAct,safe='*')
  nullify(CMO)
end subroutine MBPT2_Clean

end module MBPT2_Global
