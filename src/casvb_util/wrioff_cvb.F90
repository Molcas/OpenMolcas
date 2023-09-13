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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine wrioff_cvb(ifield,file_id,ioffset)

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
parameter(nbuf=50)
dimension ioff(nbuf)

if (ifield > nbuf) then
  write(6,*) ' ifield too large in wrioff :',ifield,nbuf
  call abend_cvb()
end if
if (tstfile_cvb(file_id)) then
  call rdi_cvb(ioff,nbuf,file_id,0)
else
  call izero(ioff,nbuf)
end if
ioff(ifield) = ioffset
call wri_cvb(ioff,nbuf,file_id,0)

return

end subroutine wrioff_cvb
