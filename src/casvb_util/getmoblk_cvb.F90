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

subroutine getmoblk_cvb(cmoblk)

use casvb_global, only: filename, nbasisq_mo, strtmo
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: cmoblk(nbasisq_mo)
integer(kind=iwp) :: iad, iadr15(15), ibf, lujob

call mkfn_cvb(strtmo,ibf)
lujob = 15
call daname_cvb(lujob,filename(ibf))
iad = 0
call iDaFile(lujob,2,iadr15,15,iad)
iad = iadr15(2)
call dDaFile(lujob,2,cmoblk,nbasisq_mo,iad)
call daclos_cvb(lujob)

return

end subroutine getmoblk_cvb
