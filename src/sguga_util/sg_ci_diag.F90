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
subroutine sg_diag_reference(SGS,CIS,EXS, &
                             nCSFs,PsiSym, &
                             TUVX_Tri,nTUVX_Tri, &
                             TU_Tri,nTU_Tri, &
                             Diag)
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: iwp, wp

   real(kind=wp), intent(out) :: Diag(nCSFs)

   real(kind=wp), allocatable :: Psi(:), Sigma(:)

   integer(kind=iwp) :: iCSF

   call mma_allocate(Psi,nCSFs,label='Psi')
   call mma_allocate(Sigma,nCSFs,label='Sigma')

   do iCSF=1,nCSFs

      Psi(:)   = Zero
      Sigma(:) = Zero

      Psi(iCSF) = One

      call sg_h_psi(SGS,CIS,EXS, &
                    Psi,nCSFs,PsiSym, &
                    Sigma, &
                    TUVX_Tri,nTUVX_Tri, &
                    TU_Tri,nTU_Tri)

      Diag(iCSF) = Sigma(iCSF)

   end do

   call mma_deallocate(Sigma)
   call mma_deallocate(Psi)

end subroutine sg_diag_reference
