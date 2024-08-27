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

subroutine Free_DeDe(Dens,TwoHam,nDens)

use k2_arrays, only: pDq, pFq, DeDe, Dq, Fq, ipOffD
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants, only: Two, Half
use stdalloc, only: mma_deallocate

implicit none
integer nDens
real*8 :: Dens(nDens), TwoHam(nDens)
integer nC, ijQ, jiQ, ij, i, j

nullify(pDq)
nullify(pFq)

if (nIrrep == 1) then
  ! symmetrize fock matrix
  ! Fix the diagonal elements of D and F
  call DScal_(nDens,Two,Dens,1)
  nc = nbas(0)
  ijq = 0
  jiq = 1-nc
  ij = 0
  do i=1,nc
    do j=1,i
      ij = ij+1
      TwoHam(ij) = Half*(Fq(ijq+j)+Fq(jiq+j*nc))
    end do
    Dens(ij) = Half*Dens(ij)
    jiq = jiq+1
    ijq = ijq+nc
  end do
  call mma_deallocate(Dq)
  call mma_deallocate(Fq)
end if

call mma_deallocate(ipOffD)
call mma_deallocate(DeDe)

return

end subroutine Free_DeDe
