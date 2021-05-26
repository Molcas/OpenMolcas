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

subroutine LDF_ComputeU(AP_QD,nD,ip_DBlocks,U)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute
!
!      U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|
!
! Note: diagonal integrals for A=B must be stored quadratically.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: AP_QD(*), nD, ip_DBlocks(nD)
real(kind=wp), intent(out) :: U(nD)
integer(kind=iwp) :: iD, AB, A, B, nA, nB, uv, ipDel, ipDB
integer(kind=iwp), external :: LDF_nBas_Atom
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

do iD=1,nD
  U(iD) = Zero
  do AB=1,NumberOfAtomPairs
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nA = LDF_nBas_Atom(A)
    nB = LDF_nBas_Atom(B)
    ipDel = AP_QD(AB)-1
    ipDB = iWork(ip_DBlocks(iD)-1+AB)-1
    do uv=1,nA*nB
      U(iD) = U(iD)+sqrt(Work(ipDel+uv))*abs(Work(ipDB+uv))
    end do
  end do
end do

end subroutine LDF_ComputeU
