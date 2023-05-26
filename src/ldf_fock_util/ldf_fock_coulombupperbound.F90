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
! Copyright (C) 2011, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombUpperBound(PrintNorm,nD,FactC,ip_DBlocks,ip_FBlocks)
!  Thomas Bondo Pedersen, January 2011.
!
!  Purpose: add the upper bound correction to the LDF Fock
!           matrix (Coulomb only):
!
!        |Fexact(uv)-F(uv)| <= sqrt[(Delta(uv)|Delta(uv))]*U
!
!        U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: PrintNorm
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: FactC(nD)
integer(kind=iwp) :: l_U, ip, l_Norm, iD, AB
real(kind=wp) :: UBFNorm
real(kind=wp), allocatable :: U(:), Norm(:)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

if (nD < 1) return
if (NumberOfAtomPairs < 1) return

l_U = nD
call mma_allocate(U,l_U,label='LDFCU')
!-tbp call LDF_GetQuadraticDiagonal(ip)
ip = ip_AP_Diag
call LDF_ComputeU(iWork(ip),nD,ip_DBlocks,U)
call LDF_Fock_CUB(iWork(ip),nD,FactC,U,ip_FBlocks)
!-tbp call LDF_FreeQuadraticDiagonal(ip)
call mma_deallocate(U)

if (PrintNorm) then
  if (NumberOfAtomPairs > 0) then
    l_Norm = NumberOfAtomPairs
    call mma_allocate(Norm,l_Norm,label='UBFNrm')
    do iD=1,nD
      call LDF_BlockMatrixNorm(iWork(ip_FBlocks(iD)),Norm)
      UBFNorm = Zero
      do AB=1,NumberOfAtomPairs
        UBFNorm = UBFNorm+Norm(AB)**2
      end do
      write(u6,'(A,I10,A,1P,D20.10,1X,A,D20.10,A)') 'Norm of Fock matrix after adding Coulomb upper bound for density',iD,':', &
                                                    sqrt(UBFNorm),'(BlockRMS=',sqrt(UBFNorm/real(NumberOfAtomPairs,kind=wp)),')'
    end do
    call xFlush(u6)
    call mma_deallocate(Norm)
  end if
end if

end subroutine LDF_Fock_CoulombUpperBound
