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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

subroutine Get_Two_Ind(IndPUVX,IndTUVX)
!***********************************************************************
! Readapted from src/fock_util/get_tuvx.f
! Return to an index in the PUVX array given
! four MO indices.

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use input_mclr, only: nAsh, nIsh, nOrb, nSym, ntAsh, ntBas
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: IndPUVX(ntBas,ntAsh,ntAsh,ntAsh), IndTUVX(ntAsh,ntAsh,ntAsh,ntAsh)
integer(kind=iwp) :: iAsh, iIsh, iIT, iIU, iIV, iIX, ijSym, iO, iOrb, iP, iPUVX, iStack, iSym, iT, iTemp, iTU, iU, iV, iVX, iX, &
                     jAsh, jO, jSym, kAsh, kl_Orb_Pairs, klSym, kO, kSym, lAsh, lMax, lO, lSym, off_Ash(nSym), off_Orb(nSym), &
                     off_PUVX(nSym)

! generate offsets

! Initialization
IndTUVX(:,:,:,:) = 0
IndPUVX(:,:,:,:) = 0

off_Orb(1) = 0
off_Ash(1) = 0
do iSym=2,nSym
  off_Orb(iSym) = off_Orb(iSym-1)+nOrb(iSym-1)
  off_Ash(iSym) = off_Ash(iSym-1)+nAsh(iSym-1)
end do

iStack = 0
do iSym=1,nSym
  off_PUVX(iSym) = iStack
  iOrb = nOrb(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = Mul(iSym,jSym)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = Mul(kSym,lSym)
        if (ijSym == klSym) then
          kl_Orb_pairs = kAsh*lAsh
          if (kSym == lSym) kl_Orb_pairs = nTri_Elem(kAsh)
          iStack = iStack+iOrb*jAsh*kl_Orb_pairs
        end if
      end do
    end do
  end do
end do

! select integrals with all 4 indices active

do iSym=1,nSym
  iOrb = nOrb(iSym)
  iAsh = nAsh(iSym)
  iIsh = nIsh(iSym)
  iPUVX = off_PUVX(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = Mul(iSym,jSym)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      lSym = Mul(ijSym,kSym)
      lAsh = nAsh(lSym)

      if ((lSym <= kSym) .and. (iAsh*jAsh*kAsh*lAsh /= 0)) then
        do iV=1,kAsh
          lMax = lAsh
          if (kSym == lSym) lMax = iV
          do iX=1,lMax
            do iU=1,jAsh
              do iP=1,iOrb
                iT = iP-iIsh
                iPUVX = iPUVX+1
                io = iP+Off_Orb(Isym)
                jo = iU+Off_Ash(Jsym)
                ko = iV+Off_Ash(ksym)
                lo = iX+Off_Ash(lsym)
                IndPUVX(io,jo,ko,lo) = iPUVX
                IndPUVX(io,jo,lo,ko) = iPUVX
                if ((iT > 0) .and. (iT <= iAsh)) then
                  iiT = iT+off_Ash(iSym)
                  iiU = iU+off_Ash(jSym)
                  if (iiU > iiT) then
                    iiT = iU+off_Ash(jSym)
                    iiU = iT+off_Ash(iSym)
                  end if

                  iTU = nTri_Elem(iiT-1)+iiU
                  iiV = iV+off_Ash(kSym)
                  iiX = iX+off_Ash(lSym)
                  if (iiX > iiV) then
                    iiV = iX+off_Ash(lSym)
                    iiX = iV+off_Ash(kSym)
                  end if
                  iVX = nTri_Elem(iiV-1)+iiX
                  if (iVX > iTU) then
                    iTemp = iTU
                    iTU = iVX
                    iVX = iTemp
                  end if
                  IndTUVX(iiT,iiU,iiV,iiX) = iPUVX
                  IndTUVX(iiT,iiU,iiX,iiV) = iPUVX
                  IndTUVX(iiU,iiT,iiV,iiX) = iPUVX
                  IndTUVX(iiU,iiT,iiX,iiV) = iPUVX
                  IndTUVX(iiV,iiX,iiT,iiU) = iPUVX
                  IndTUVX(iiV,iiX,iiU,iiT) = iPUVX
                  IndTUVX(iiX,iiV,iiT,iiU) = iPUVX
                  IndTUVX(iiX,iiV,iiU,iiT) = iPUVX
                end if
              end do
            end do
          end do
        end do
      end if
    end do
  end do
end do

end subroutine Get_Two_Ind
