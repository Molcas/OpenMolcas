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

implicit none
logical PrintNorm
integer nD
real*8 FactC(nD)
integer ip_DBlocks(nD)
real*8 UBFNorm(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer LDF_nBas_Atom
external LDF_nBas_Atom

integer ip
integer ip_U, l_U
integer iD
integer AB
integer nAB
integer ipDel
integer uv

integer i, j
integer AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

if (nD < 1) return
if (NumberOfAtomPairs < 1) return

!-tbp Call LDF_GetQuadraticDiagonal(ip)
ip = ip_AP_Diag
l_U = nD
call GetMem('CUBNrmU','Allo','Real',ip_U,l_U)
call LDF_ComputeU(ip,nD,ip_DBlocks,Work(ip_U))
do iD=1,nD
  UBFNorm(iD) = 0.0d0
  do AB=1,NumberOfAtomPairs
    nAB = LDF_nBas_Atom(AP_Atoms(1,AB))*LDF_nBas_Atom(AP_Atoms(2,AB))
    ipDel = iWork(ip-1+AB)-1
    do uv=1,nAB
      UBFNorm(iD) = UBFNorm(iD)+Work(ipDel+uv)
    end do
  end do
  UBFNorm(iD) = FactC(iD)*Work(ip_U-1+iD)*sqrt(UBFNorm(iD))
end do
call GetMem('CUBNrmU','Free','Real',ip_U,l_U)
!-tbp Call LDF_FreeQuadraticDiagonal(ip)

if (PrintNorm) then
  do iD=1,nD
    write(6,'(A,I10,A,1P,D20.10,1X,A,D20.10,A)') 'Norm of upper bound Coulomb error for density',iD,':',UBFNorm(iD),'(BlockRMS=', &
                                                 sqrt(UBFNorm(iD)**2/dble(NumberOfAtomPairs)),')'
  end do
  call xFlush(6)
end if

end subroutine LDF_Fock_CoulombUpperBoundNorm
