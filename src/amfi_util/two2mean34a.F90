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

subroutine two2mean34a(carteSO,carteOO,occup,AOcoeffs,onecart,ncontmf,norbsum,noccorb,sameorb)

use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
#include "para.fh"
integer(kind=iwp) :: ncontmf, norbsum, noccorb
real(kind=wp) :: carteSO(norbsum,ncontmf,norbsum,ncontmf), carteOO(norbsum,ncontmf,norbsum,ncontmf), occup(*), &
                 AOcoeffs(MxcontL,*), onecart(MxcontL,MxcontL)
logical(kind=iwp) :: sameorb
integer(kind=iwp) :: icartleft, icartright, irun, jrun, Mrun
real(kind=wp) :: coeff

if (sameorb) then
  do icartleft=1,norbsum
    do icartright=1,norbsum
      coeff = Zero
      do Mrun=1,noccorb
        coeff = coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*AOcoeffs(icartright,Mrun)
      end do
      coeff = Half*coeff
      do irun=1,ncontmf
        do jrun=1,ncontmf
          onecart(irun,jrun) = onecart(irun,jrun)+coeff*carteSO(icartleft,irun,icartright,jrun)
        end do
      end do
    end do
  end do
else
  do icartleft=1,norbsum
    do icartright=1,norbsum
      coeff = Zero
      do Mrun=1,noccorb
        coeff = coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*AOcoeffs(icartright,Mrun)
      end do
      coeff = Half*coeff
      do irun=1,ncontmf
        do jrun=1,ncontmf
          onecart(irun,jrun) = onecart(irun,jrun)+ &
                               coeff*(carteSO(icartleft,irun,icartright,jrun)+Two*carteOO(icartleft,irun,icartright,jrun))
        end do
      end do
    end do
  end do
end if

return

end subroutine two2mean34a
