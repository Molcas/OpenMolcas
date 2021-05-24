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
! Copyright (C) 2014, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombOnly_XIDI(Mode,tau,nD,FactC,ip_DBlocks,ip_FBlocks)
! Thomas Bondo Pedersen, June 2014.
!
! Compute correction from errors in the integral diagonal blocks,
! adding them to the blocked Fock matrix.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Mode, nD, ip_DBlocks(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: tau, FactC(nD)
logical(kind=iwp) :: IPI_set_here
integer(kind=iwp) :: TaskListID, AB, A, B, nAB, ip_Int, l_Int, n_Int, iD, ipD, ipF
logical(kind=iwp), external :: Rsv_Tsk, LDF_IntegralPrescreeningInfoIsSet
integer(kind=iwp), external :: LDF_nBas_Atom
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
! statement function
integer(kind=iwp) :: i, j, AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

if (.not. LDF_IntegralPrescreeningInfoIsSet()) then
  call LDF_SetIntegralPrescreeningInfo()
  IPI_set_here = .true.
else
  IPI_set_here = .false.
end if

call Init_Tsk(TaskListID,NumberOfAtomPairs)
do while (Rsv_Tsk(TaskListID,AB))
  A = AP_Atoms(1,AB)
  B = AP_Atoms(2,AB)
  nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
  if (nAB > 0) then
    n_Int = nAB**2
    l_Int = 2*n_Int
    call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
    call LDF_ComputeValenceIntegrals(AB,AB,n_Int,Work(ip_Int))
    call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,AB,n_Int,Work(ip_Int+n_Int))
    call dAXPY_(n_Int,-One,Work(ip_Int+n_Int),1,Work(ip_Int),1)
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+AB)
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      call dGeMV_('N',nAB,nAB,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,One,Work(ipF),1)
    end do
    call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
  end if
end do
call Free_Tsk(TaskListID)

if (IPI_set_here) then
  call LDF_UnsetIntegralPrescreeningInfo()
end if

end subroutine LDF_Fock_CoulombOnly_XIDI
