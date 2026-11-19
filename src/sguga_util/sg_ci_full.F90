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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

subroutine sg_ci_full(SGS,CIS,EXS,nCSFs,CSFSym, &
                      TUVX_Tri,nTUVX_Tri, &
                      TU_Tri,nTU_Tri, &
                      Full)

  use sguga, only: SGStruct, CIStruct, EXStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp

  implicit none

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in) :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Full(nCSFs,nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Sigma(:)

  integer(kind=iwp) :: iCSF, jCSF

  call mma_allocate(Sigma,nCSFs,Label='Sigma')
  call mma_allocate(Psi,nCSFs,Label='Psi')

  do iCSF=1,nCSFs

     Psi(:)   = Zero
     Sigma(:) = Zero

     Psi(iCSF) = One

     ! Construct H|I> = Sum_J <J|H|I> |J>

     call sg_h_psi(SGS,CIS,EXS, &
                   Psi,nCSFs,CSFSym, &
                   Sigma, &
                   TUVX_Tri,nTUVX_Tri, &
                   TU_Tri,nTU_Tri)

     Do jCSF=1,nCSFs
        Full(jCSF,iCSF) = Sigma(jCSF)
     End Do

  end do
  call RecPrt('CI matrix in CSF basis',' ',Full,nCSFs,nCSFs)

  call mma_deallocate(Sigma)
  call mma_deallocate(Psi)

end subroutine sg_ci_full
