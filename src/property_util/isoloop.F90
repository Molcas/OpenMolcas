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

subroutine isoloop(double)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
logical(kind=iwp), intent(in) :: double
integer(kind=iwp) :: n, nAtoms_All, nTemp
character(len=2) :: Element(MxAtom)
real(kind=wp), allocatable :: Coor(:,:), T(:)

!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,*)
call CollapseOutput(1,'   Isotopic shifts:')
write(u6,'(3X,A)') '   ----------------'
write(u6,*)

call Get_nAtoms_All(nAtoms_All)
call mma_allocate(Coor,3,nAtoms_All)
call Get_Coord_All(Coor,nAtoms_All)
call Get_Name_All(Element)

write(u6,*)
write(u6,*)
write(u6,*) '****************************************'
write(u6,*) '*                                      *'
write(u6,*) '* Isotope shifted frequencies in cm-1  *'
write(u6,*) '*                                      *'
write(u6,*) '****************************************'
write(u6,*)
n = 3*nAtoms_All
nTemp = 6*2*n**2
call mma_allocate(T,nTemp,label='ISOLOOP')
call Isotop_i(n,Element,nAtoms_All,T,nTemp,Coor,double)
call mma_deallocate(T)
call mma_deallocate(Coor)

call CollapseOutput(0,'   Isotopic shifts:')
write(u6,*)

return

end subroutine isoloop
