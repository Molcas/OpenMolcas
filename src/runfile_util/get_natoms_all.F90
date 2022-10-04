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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Get_nAtoms_All
!
!> @author R. Lindh
!>
!> @details
!> Get number of all atoms (not only symmetry unique) from RUNFILE.
!>
!> @param[out] nAtoms_All Number of all atoms in the molecule
!***********************************************************************

subroutine Get_nAtoms_All(nAtoms_All)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nAtoms_All
#include "WrkSpc.fh"
integer(kind=iwp) :: ipCoord, nAtoms

call Get_iScalar('Unique atoms',nAtoms)
call Allocate_Work(ipCoord,3*nAtoms)
call Get_dArray('Unique Coordinates',Work(ipCoord),3*nAtoms)
call Get_nAtoms_All_(Work(ipCoord),nAtoms,nAtoms_all)
call Free_Work(ipCoord)

!write(u6,*) 'nAtoms_All=',nAtoms_All
!write(u6,*) 'Exit Get_nAtoms_All'
return

end subroutine Get_nAtoms_All
