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

subroutine OneBas(Label)
!***********************************************************************
!                                                                      *
!     Change nBas in OneDat.fh                                         *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
integer IntBas(8)
character*(*) Label
#include "OneDat.fh"

if (Label == 'CONT') then
  call Get_iArray('nBas',IntBas,nSym)
else if (Label == 'PRIM') then
  call Get_iArray('nBas_Prim',IntBas,nSym)
else
  write(6,*) 'OneBas: Illegal Label value!'
  write(6,*) 'Value: ',Label
  call Abend()
end if
call ICopy(nSym,IntBas,1,nBas,1)

return

end subroutine OneBas
