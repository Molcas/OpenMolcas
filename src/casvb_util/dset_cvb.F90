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

subroutine dset_cvb(iorbrel,ifxorb,ifxstr,idelstr,iorts,irots,izeta)

use casvb_global, only: mxorb_cvb, ndimrel, ndrot, nfxvb, norb, norbrel, nort, nsyme, nzrvb, plc_const, recinp, sym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iorbrel(ndimrel), ifxorb(mxorb_cvb), ifxstr(nfxvb), idelstr(nzrvb), iorts(*), irots(*), izeta(*)
integer(kind=iwp) :: ioffs

! Check if any molecular interaction constraints:
!if (ploc) call plcconst_plc()
sym = ((norbrel > 0) .or. (nort > 0) .or. plc_const)
call rdioff_cvb(9,recinp,ioffs)
call wris_cvb(iorbrel,ndimrel,recinp,ioffs)
call wrioff_cvb(10,recinp,ioffs)
call wrioff_cvb(11,recinp,ioffs)
call wris_cvb(ifxorb,norb,recinp,ioffs)
call wrioff_cvb(12,recinp,ioffs)
call wris_cvb(ifxstr,nfxvb,recinp,ioffs)
call wrioff_cvb(13,recinp,ioffs)
call wris_cvb(idelstr,nzrvb,recinp,ioffs)
call wrioff_cvb(14,recinp,ioffs)
call wris_cvb(iorts,2*nort,recinp,ioffs)
call wrioff_cvb(15,recinp,ioffs)
call wris_cvb(irots,2*ndrot,recinp,ioffs)
call wrioff_cvb(16,recinp,ioffs)
call wris_cvb(izeta,nsyme,recinp,ioffs)
call wrioff_cvb(17,recinp,ioffs)

return

end subroutine dset_cvb
