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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_FCI(UsePartPermSym,nD,FactC,ip_DBlocks,ip_FBlocks)
! Thomas Bondo Pedersen, October 2010.
!
! Purpose: Compute Coulomb contribution to Fock matrix using
!          conventional integrals (debug code).

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: UsePartPermSym
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: FactC(nD)
integer(kind=iwp) :: AB, CD, A, B, C, D, nAB, nCD, ip_Int, l_Int, iD, ipD, ipF
integer(kind=iwp), external :: LDF_nBas_Atom
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
!statement function
integer(kind=iwp) :: i, j, AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

if (UsePartPermSym) then ! use particle permutation symmetry
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,AB-1
      C = AP_Atoms(1,CD)
      D = AP_Atoms(2,CD)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
      call LDF_ComputeValenceIntegrals(AB,CD,l_Int,Work(ip_Int))
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,One,Work(ipF),1)
      end do
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+AB)
        ipF = iWork(ip_FBlocks(iD)-1+CD)
        call dGeMV_('T',nAB,nCD,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,One,Work(ipF),1)
      end do
      call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
    end do
    l_Int = nAB**2
    call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
    call LDF_ComputeValenceIntegrals(AB,AB,l_Int,Work(ip_Int))
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+AB)
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      call dGeMV_('N',nAB,nAB,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,One,Work(ipF),1)
    end do
    call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
  end do
else
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,NumberOfAtomPairs
      C = AP_Atoms(1,CD)
      D = AP_Atoms(2,CD)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
      call LDF_ComputeValenceIntegrals(AB,CD,l_Int,Work(ip_Int))
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),Work(ip_Int),nAB,Work(ipD),1,One,Work(ipF),1)
      end do
      call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
    end do
  end do
end if

end subroutine LDF_FCI
