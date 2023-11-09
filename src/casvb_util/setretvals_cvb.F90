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

subroutine setretvals_cvb(esym,n_iter)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n_iter
real(kind=wp), intent(in) :: esym(*)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

if (nac == 0) then
  ener(1,iter) = emy
else
  ener(1:lroots,iter) = esym(stsym)
end if
iterci = n_iter

return

end subroutine setretvals_cvb
