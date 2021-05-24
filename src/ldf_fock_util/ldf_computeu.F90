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

subroutine LDF_ComputeU(ip_AP_QD,nD,ip_DBlocks,U)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute
!
!      U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|
!
! Note: diagonal integrals for A=B must be stored quadratically.

implicit none
integer ip_AP_QD
integer nD
integer ip_DBlocks(nD)
real*8 U(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer LDF_nBas_Atom
external LDF_nBas_Atom

integer iD
integer AB
integer A, B
integer nA, nB
integer uv
integer ipDel, ipDB

integer i, j
integer ip_Delta
integer AP_Atoms
ip_Delta(i) = iWork(ip_AP_QD-1+i)
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

do iD=1,nD
  U(iD) = 0.0d0
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nA = LDF_nBas_Atom(A)
    nB = LDF_nBas_Atom(B)
    ipDel = ip_Delta(AB)-1
    ipDB = iWork(ip_DBlocks(iD)-1+AB)-1
    do uv=1,nA*nB
      U(iD) = U(iD)+sqrt(Work(ipDel+uv))*abs(Work(ipDB+uv))
    end do
  end do
end do

end subroutine LDF_ComputeU
