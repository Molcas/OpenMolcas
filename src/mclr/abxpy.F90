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

use MCLR_Data, only: ipMO, NA
use input_mclr, only: nSym, nAsh, nIsh, nOrb

implicit none
integer idsym
real*8 Array1(*), Array2(*)
integer i, j, itri
integer iS, jS, kS, lS, ijS
integer iA, jA, kA, lA
integer iAsh, jAsh, kAsh, lAsh
integer iiA, ij, kl, ijkl, ip1
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

do iS=1,nSym
  do iA=1,Nash(is)
    iAsh = nA(is)+iA
    iiA = nIsh(is)+iA
    do jS=1,nSym
      ijs = ieor(is-1,js-1)+1
      do jA=1,Nash(js)
        jAsh = nA(js)+jA
        ij = itri(iash,jash)
        if (iAsh >= jash) then
          do kS=1,nSym
            do kA=1,Nash(ks)
              kAsh = nA(ks)+kA
              ls = ieor(ieor(kS-1,ijs-1),idsym-1)+1
              do lA=1,Nash(ls)
                lAsh = nA(ls)+lA
                if (kAsh >= lash) then
                  kl = itri(kAsh,lash)
                  if (ij >= kl) then
                    ijkl = itri(ij,kl)
                    ip1 = ipMO(js,ks,ls)+nOrb(is)*nAsh(js)*((lA-1)*nAsh(kS)+kA-1)+nOrb(is)*(ja-1)+iia-1
                    Array2(ijkl) = array1(ip1)+array2(ijkl)
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
