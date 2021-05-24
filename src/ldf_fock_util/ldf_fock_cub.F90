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

subroutine LDF_Fock_CUB(ip_AP_QD,nD,FactC,U,ip_FBlocks)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute
!
!      sqrt[(Delta(uv)|Delta(uv))]*U
!
!      and add to Fock matrix.
!
! Note: diagonal integrals for A=B must be stored quadratically.

implicit none
integer ip_AP_QD
integer nD
real*8 FactC(nD)
real*8 U(nD)
integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer LDF_nBas_Atom
external LDF_nBas_Atom

integer iD
integer AB
integer A, B
integer nA, nB
integer uv
integer ipDel, ipFB

real*8 UU

integer i, j
integer ip_Delta
integer AP_Atoms
ip_Delta(i) = iWork(ip_AP_QD-1+i)
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

do iD=1,nD
  UU = FactC(iD)*U(iD)
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nA = LDF_nBas_Atom(A)
    nB = LDF_nBas_Atom(B)
    ipDel = ip_Delta(AB)-1
    ipFB = iWork(ip_FBlocks(iD)-1+AB)-1
    do uv=1,nA*nB
      Work(ipFB+uv) = Work(ipFB+uv)+sqrt(Work(ipDel+uv))*UU
    end do
  end do
end do

end subroutine LDF_Fock_CUB
