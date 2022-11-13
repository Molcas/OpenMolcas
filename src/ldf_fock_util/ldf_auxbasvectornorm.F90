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

subroutine LDF_AuxBasVectorNorm(V,Norm)
! Thomas Bondo Pedersen, December 2010.
!
! Purpose: compute Frobenius norm of aux bas vector.

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: V(*)
real(kind=wp), intent(_OUT_) :: Norm(*)
integer(kind=iwp) :: iAtom, nAtom, iAtomPair, ip, l
integer(kind=iwp), external :: LDF_nAtom, LDF_nBasAux_Atom
real(kind=wp), external :: ddot_
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

nAtom = LDF_nAtom()
do iAtom=1,nAtom
  l = LDF_nBasAux_Atom(iAtom)
  ip = V(iAtom)
  Norm(iAtom) = sqrt(dDot_(l,Work(ip),1,Work(ip),1))
end do
do iAtomPair=1,NumberOfAtomPairs
  l = iWork(ip_AP_2CFunctions-1+2*(iAtomPair-1)+1)
  ip = V(nAtom+iAtomPair)
  Norm(nAtom+iAtomPair) = sqrt(dDot_(l,Work(ip),1,Work(ip),1))
end do

end subroutine LDF_AuxBasVectorNorm
