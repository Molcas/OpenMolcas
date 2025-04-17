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

subroutine Pickmo_MCLR(rmo,rmoaa,idsym)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMO, nA
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rmo(*)
real(kind=wp), intent(_OUT_) :: rmoaa(*)
integer(kind=iwp), intent(in) :: idsym
integer(kind=iwp) :: iA, iAA, ijAA, ijkl, ipi, iS, jA, jAA, jS, kA, kAA, klAA, kS, lA, lAA, lS

do iS=1,nSym
  do jS=1,iS
    do kS=1,is
      ls = Mul(Mul(is,js),Mul(ks,idsym))
      if (ls <= ks) then
        do iA=1,nAsh(is)
          iAA = iA+nA(is)
          do jA=1,nAsh(js)
            jAA = jA+nA(js)
            ijAA = iTri(iAA,jAA)
            do kA=1,nAsh(ks)
              kAA = kA+nA(ks)
              do lA=1,nAsh(ls)
                lAA = lA+nA(ls)
                klAA = iTri(kAA,lAA)
                if (ijAA >= klAA) then
                  ijkl = iTri(ijAA,klAA)
                  ipi = ipMO(js,ks,ls)+nIsh(is)+iA-1+nOrb(is)*(jA-1)+nOrb(is)*nAsh(js)*(kA-1)+nOrb(is)*nAsh(js)*nAsh(ks)*(lA-1)
                  !ipi = ipMO(js,ks,ls)+(nOrb(is)*(jA-1+nAsh(js)*(kA-1+nAsh(ks)*(lA-1))))+nish(is)+iA-1
                  rmoaa(ijkl) = rmo(ipi)
                end if
              end do
            end do
          end do
        end do
      end if
    end do
  end do
end do

end subroutine Pickmo_MCLR
