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

subroutine Get_Two_Ind(Ind_PUVX,IndTUVX)
!***********************************************************************
! Readapted from src/fock_util/get_tuvx.f
! Return to an index in the PUVX array given
! four MO indices.

use input_mclr, only: nSym, ntBas, ntAsh, nAsh, nIsh, nOrb

implicit none
! Output
integer, dimension(ntBas,ntAsh,ntAsh,ntAsh) :: Ind_PUVX
integer, dimension(ntAsh,ntAsh,ntAsh,ntAsh) :: IndTUVX
! Auxiliaries
integer, dimension(nSym) :: off_Ash, off_PUVX, off_Orb
integer lOrb, kOrb, jOrb, iOrb, iStack, iSym, jSym, jAsh, ijSym, kSym, kAsh, lSym, lAsh, klSym, kl_Orb_Pairs, iAsh, iIsh, iPUVX, &
        iV, lMax, iX, iU, iP, iT, iO, jO, kO, lO, iIT, iIU, iTU, iIV, iIX, iVX, iTemp
! Statement function
integer i, iTri
iTri(i) = (i*i-i)/2

! generate offsets

! Initialization
do lOrb=1,ntAsh
  do KOrb=1,ntAsh
    do JOrb=1,ntAsh
      do iOrb=1,ntAsh
        IndTUVX(iOrb,jOrb,kOrb,lOrb) = 0
        Ind_PUVX(iOrb,jOrb,kOrb,lOrb) = 0
      end do
      do iOrb=ntAsh+1,ntBas
        Ind_PUVX(iOrb,jOrb,kOrb,lOrb) = 0
      end do
    end do
  end do
end do

iStack = 0
do iSym=1,nSym
  off_Orb(iSym) = iStack
  iStack = iStack+nOrb(iSym)
end do

iStack = 0
do iSym=1,nSym
  off_Ash(iSym) = iStack
  iStack = iStack+nAsh(iSym)
end do

iStack = 0
do iSym=1,nSym
  off_PUVX(iSym) = iStack
  iOrb = nOrb(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = 1+ieor(iSym-1,jSym-1)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      do lSym=1,kSym
        lAsh = nAsh(lSym)
        klSym = 1+ieor(kSym-1,lSym-1)
        if (ijSym == klSym) then
          kl_Orb_pairs = kAsh*lAsh
          if (kSym == lSym) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
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
    ijSym = 1+ieor(iSym-1,jSym-1)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      lSym = 1+ieor(ijSym-1,kSym-1)
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
                Ind_PUVX(io,jo,ko,lo) = iPUVX
                Ind_PUVX(io,jo,lo,ko) = iPUVX
                if ((iT > 0) .and. (iT <= iAsh)) then
                  iiT = iT+off_Ash(iSym)
                  iiU = iU+off_Ash(jSym)
                  if (iiU > iiT) then
                    iiT = iU+off_Ash(jSym)
                    iiU = iT+off_Ash(iSym)
                  end if

                  iTU = iiU+iTri(iiT)
                  iiV = iV+off_Ash(kSym)
                  iiX = iX+off_Ash(lSym)
                  if (iiX > iiV) then
                    iiV = iX+off_Ash(lSym)
                    iiX = iV+off_Ash(kSym)
                  end if
                  iVX = iiX+iTri(iiV)
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
