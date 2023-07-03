!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

integer function nFgBond(iAtom,nAtoms,nMax,iTabBonds,nBondMax,iTabAtoms)

implicit real*8(a-h,o-z)
integer iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
nFgBond = 0
nn = iTabAtoms(1,0,iAtom)
do i=1,nn
  if (iTabBonds(3,iTabAtoms(2,i,iAtom)) == Fragments_Bond) nFgBond = nFgBond+1
end do

return

end function nFgBond
