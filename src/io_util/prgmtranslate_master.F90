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
! Copyright (C) 2001-2016, Valera Veryazov                             *
!***********************************************************************

subroutine prgmtranslate_master(in,out,lout)

character*(*) in, out
integer Strnln
external Strnln

lin = Strnln(in)
out = ' '
if (index(in,'/') /= 0) then
  ! just in case if we processing translated name!
  out = in
  lout = lin
else
  call prgmtranslatec(in,lin,out,lout,0)
end if
out = out(1:lout)

return

end subroutine prgmtranslate_master
