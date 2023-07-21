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

logical function Cho_Rsv_Tsk(ID,i)

implicit none
integer ID, i
#include "cho_para_info.fh"
logical Rsv_Tsk
external Rsv_Tsk
logical Cho_Rsv_Tsk_
external Cho_Rsv_Tsk_

if (Cho_Real_Par) then
  Cho_Rsv_Tsk = Rsv_Tsk(ID,i)
else
  Cho_Rsv_Tsk = Cho_Rsv_Tsk_(ID,i)
end if

end function Cho_Rsv_Tsk
