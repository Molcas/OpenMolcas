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

subroutine Print_Input(Title,nSym,nBas,wSet,nSet)

implicit real*8(a-h,o-z)

#include "mxdm.fh"
#include "mxave.fh"

dimension wSet(nSet), nBas(MxSym)
character*72 Title
parameter(lPaper=132)

call Banner(Title,1,lPaper-7)
write(6,*)
write(6,*)
write(6,*) 'Number of symmetries:',nSym
write(6,*) 'Basis functions:',(nBas(iSym),iSym=1,nSym)
write(6,*) 'Normalized weights:',(wSet(iS),iS=1,nSet)

return

end subroutine Print_Input
