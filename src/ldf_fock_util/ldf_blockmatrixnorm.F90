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

subroutine LDF_BlockMatrixNorm(Blocks,Norm)
! Thomas Bondo Pedersen, December 2010.
!
! Purpose: compute Frobenius norm of block matrix.

use Definitions, only: wp, iwp

implicit none
#include "ldf_atom_pair_info.fh"
integer(kind=iwp), intent(in) :: Blocks(*)
real(kind=wp), intent(out) :: Norm(NumberOfAtomPairs)
integer(kind=iwp) :: iAtomPair, iAtom, jAtom, ip, l
integer(kind=iwp), external :: LDF_nBas_Atom
real(kind=wp), external :: ddot_
#include "WrkSpc.fh"

do iAtomPair=1,NumberOfAtomPairs
  iAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+1)
  jAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+2)
  l = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  ip = Blocks(iAtomPair)
  Norm(iAtomPair) = sqrt(dDot_(l,Work(ip),1,Work(ip),1))
end do

end subroutine LDF_BlockMatrixNorm
