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

subroutine biksmain_cvb(aikcof,bikcof,nel,nalf,ndet,ifns,kbasis,share,iprint)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nel, nalf, ndet, ifns, kbasis, iprint
real(kind=wp), intent(inout) :: aikcof(ndet,ifns)
real(kind=wp), intent(out) :: bikcof(ndet,ifns)
logical(kind=iwp) :: share
integer(kind=iwp) :: nbet, nswpdim

if ((nel == 0) .and. (kbasis /= 6)) then
  bikcof(1,1) = One
  aikcof(1,1) = One
  return
end if

nbet = nel-nalf
nswpdim = 2**nbet

call rumer_cvb(bikcof,nel,nalf,nbet,ndet,ifns,kbasis,iprint,nswpdim)

if ((kbasis == 1) .or. (kbasis == 5)) call kotani_cvb(bikcof,ndet,ifns)

if (kbasis == 5) call projspn_cvb(bikcof,nel,nalf,nbet,ndet,ifns)

if (kbasis == 2) call serber_cvb(bikcof,nel,nalf,nbet,ndet,ifns)

call aikcof_cvb(aikcof,bikcof,ndet,ifns,kbasis,share)

return

end subroutine biksmain_cvb
