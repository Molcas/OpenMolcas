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

subroutine GetNextFfield(nextfld,fldname,nOrdOpf,ncmp,force,lforce)

implicit real*8(A-H,O-Z)
character*(*) fldname
dimension force(lforce)

nextfld = 0
return
! Avoid unused argument warnings
if (.false.) then
  call Unused_character(fldname)
  call Unused_integer(nOrdOpf)
  call Unused_integer(ncmp)
  call Unused_real_array(force)
end if

end subroutine GetNextFfield
