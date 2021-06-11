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

subroutine LDF_Fock_CoulombUpperBoundNorm(PrintNorm,nD,FactC,ip_DBlocks,UBFNorm)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute norm of the upper bound to the Coulomb Fock
!          matrix error in LDF.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: PrintNorm
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD)
real(kind=wp), intent(in) :: FactC(nD)
real(kind=wp), intent(out) :: UBFNorm(nD)
integer(kind=iwp) :: ip, l_U, iD, AB, A, B, nAB, ipDel, uv
real(kind=wp), allocatable :: U(:)
integer(kind=iwp), external :: LDF_nBas_Atom
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

if (nD < 1) return
if (NumberOfAtomPairs < 1) return

!-tbp Call LDF_GetQuadraticDiagonal(ip)
ip = ip_AP_Diag
l_U = nD
call mma_allocate(U,l_U,label='CUBNrmU')
call LDF_ComputeU(iWork(ip),nD,ip_DBlocks,U)
do iD=1,nD
  UBFNorm(iD) = Zero
  do AB=1,NumberOfAtomPairs
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    ipDel = iWork(ip-1+AB)-1
    do uv=1,nAB
      UBFNorm(iD) = UBFNorm(iD)+Work(ipDel+uv)
    end do
  end do
  UBFNorm(iD) = FactC(iD)*U(iD)*sqrt(UBFNorm(iD))
end do
call mma_deallocate(U)
!-tbp Call LDF_FreeQuadraticDiagonal(ip)

if (PrintNorm) then
  do iD=1,nD
    write(u6,'(A,I10,A,1P,D20.10,1X,A,D20.10,A)') 'Norm of upper bound Coulomb error for density',iD,':',UBFNorm(iD),'(BlockRMS=', &
                                                  sqrt(UBFNorm(iD)**2/real(NumberOfAtomPairs,kind=wp)),')'
  end do
  call xFlush(u6)
end if

end subroutine LDF_Fock_CoulombUpperBoundNorm
