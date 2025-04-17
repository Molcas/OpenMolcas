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

subroutine ABXpY(Array1,Array2,idsym)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMO, NA
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Array1(*)
real(kind=wp), intent(inout) :: Array2(*)
integer(kind=iwp), intent(in) :: idsym
integer(kind=iwp) :: iA, iAsh, iiA, ij, ijkl, ijS, ip1, iS, jA, jAsh, jS, kA, kAsh, kl, kS, lA, lAsh, lS

do iS=1,nSym
  do iA=1,Nash(is)
    iAsh = nA(is)+iA
    iiA = nIsh(is)+iA
    do jS=1,nSym
      ijs = Mul(is,js)
      do jA=1,Nash(js)
        jAsh = nA(js)+jA
        ij = iTri(iash,jash)
        if (iAsh >= jash) then
          do kS=1,nSym
            do kA=1,Nash(ks)
              kAsh = nA(ks)+kA
              ls = Mul(Mul(kS,ijs),idsym)
              do lA=1,Nash(ls)
                lAsh = nA(ls)+lA
                if (kAsh >= lash) then
                  kl = iTri(kAsh,lash)
                  if (ij >= kl) then
                    ijkl = iTri(ij,kl)
                    ip1 = ipMO(js,ks,ls)+nOrb(is)*nAsh(js)*((lA-1)*nAsh(kS)+kA-1)+nOrb(is)*(ja-1)+iia-1
                    Array2(ijkl) = Array2(ijkl)+Array1(ip1)
                  end if
                end if
              end do
            end do
          end do
        end if
      end do
    end do
  end do
end do

end subroutine ABXpY
