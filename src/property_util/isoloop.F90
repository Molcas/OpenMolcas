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

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
character*2 Element(MxAtom)
logical double

!                                                                      *
!***********************************************************************
!                                                                      *
write(6,*)
call CollapseOutput(1,'   Isotopic shifts:')
write(6,'(3X,A)') '   ----------------'
write(6,*)

call Get_nAtoms_All(nAtoms_All)
call Allocate_Work(ipCoor,3*nAtoms_All)
call Get_Coord_All(Work(ipCoor),nAtoms_All)
call Get_Name_All(Element)

write(6,*)
write(6,*)
write(6,*) '****************************************'
write(6,*) '*                                      *'
write(6,*) '* Isotope shifted frequencies in cm-1  *'
write(6,*) '*                                      *'
write(6,*) '****************************************'
write(6,*)
n = 3*nAtoms_All
nTemp = 6*2*n**2
call Getmem('ISOLOOP','ALLO','REAL',ipT,nTemp)
call Isotop_i(n,Element,nAtoms_All,Work(ipT),nTemp,Work(ipCoor),double)
call Getmem('ISOLOOP','FREE','REAL',ipT,nTemp)
call Free_Work(ipCoor)

call CollapseOutput(0,'   Isotopic shifts:')
write(6,*)

return

end subroutine isoloop
