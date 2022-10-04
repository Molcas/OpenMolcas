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

subroutine Get_nadc(Grad1,Grad2,NADC,nGrad)

implicit real*8(a-h,o-z)
#include "SysDef.fh"

real*8 Grad1(nGrad), Grad2(nGrad), NADC(nGrad)
character*24 Label
Label = 'Grad State1'
call Get_dArray(Label,Grad1,nGrad)
Label = 'Grad State2'
call Get_dArray(Label,Grad2,nGrad)
Label = 'NADC'
call Get_dArray(Label,NADC,nGrad)

return

end subroutine Get_nadc
