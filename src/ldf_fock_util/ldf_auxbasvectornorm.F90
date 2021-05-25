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

subroutine LDF_AuxBasVectorNorm(ip_V,Norm)
! Thomas Bondo Pedersen, December 2010.
!
! Purpose: compute Frobenius norm of aux bas vector.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ip_V
real(kind=wp), intent(out) :: Norm(*)
integer iAtom, nAtom, iAtomPair, ip0, ip, l
integer(kind=iwp), external :: LDF_nAtom, LDF_nBasAux_Atom
real(kind=wp), external :: ddot_
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
! statement function
integer i, j, AP_2CFunctions
AP_2CFunctions(i,j) = iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

nAtom = LDF_nAtom()
ip0 = ip_V-1
do iAtom=1,nAtom
  l = LDF_nBasAux_Atom(iAtom)
  ip = iWork(ip0+iAtom)
  Norm(iAtom) = sqrt(dDot_(l,Work(ip),1,Work(ip),1))
end do
ip0 = ip0+nAtom
do iAtomPair=1,NumberOfAtomPairs
  l = AP_2CFunctions(1,iAtomPair)
  ip = iWork(ip0+iAtomPair)
  Norm(nAtom+iAtomPair) = sqrt(dDot_(l,Work(ip),1,Work(ip),1))
end do

end subroutine LDF_AuxBasVectorNorm
