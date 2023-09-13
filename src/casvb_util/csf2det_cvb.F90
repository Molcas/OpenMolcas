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

subroutine csf2det_cvb(vec,detvec,isym_loc,iWay)

use csfbas, only: cts, kdtoc
implicit real*8(a-h,o-z)
#include "ciinfo.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
dimension vec(*), detvec(*)

if (iWay == 1) then
  if (nac == 0) then
    detvec(1) = vec(1)
    return
  end if

  jCopy = 0
  call csdtvc(vec,detvec,iway,work(kdtoc),cts,isym_loc,jcopy)
else if (iWay == 2) then
  if (nac == 0) then
    vec(1) = detvec(1)
    return
  end if

  jCopy = 0
  call csdtvc(vec,detvec,iway,work(kdtoc),cts,isym_loc,jcopy)
end if

return

end subroutine csf2det_cvb
