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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Get_TUVX(PUVX,TUVX)
!***********************************************************************
!                                                                      *
!     extract the subset of ERIs (tu!vx) from the set (pu!vx)          *
!                                                                      *
!     calling arguments:                                               *
!     PUVX    : input, array of real                                   *
!               ERIs with indices (pu!vx)                              *
!               (three active, one general)                            *
!     TUVX    : output, array of real                                  *
!               ERIs with indices (tu!vx)                              *
!               (four active)                                          *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: PUVX(*)
real(kind=wp), intent(_OUT_) :: TUVX(*)
#include "rasdim.fh"
#include "general.fh"
integer(kind=iwp) :: iAsh, iIsh, iiT, iiU, iiV, iiX, ijSym, iOrb, iP, iPUVX, iStack, iSym, iT, iTemp, iTU, iTUVX, iU, iV, iVX, iX, &
                     jAsh, jSym, kAsh, kl_Orb_pairs, klSym, kSym, lAsh, lMax, lSym, off_Ash(mxSym), off_PUVX(mxSym)

! generate offsets
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

        !write(LF,*) 'sym(p,w,x,y),offset= ',isym,jsym,ksym,lsym,iPUVX
        !call recprt('(pw|xy)','(5ES16.8)',PUVX(iPUVX+1),iorb*jAsh,kAsh*lAsh+Min(ijSym-2,0)*kAsh*(lAsh-1)/2)

        do iV=1,kAsh
          lMax = lAsh
          if (kSym == lSym) lMax = iV
          do iX=1,lMax
            do iU=1,jAsh
              do iP=1,iOrb
                iT = iP-iIsh
                iPUVX = iPUVX+1
                if ((iT > 0) .and. (iT <= iAsh)) then
                  iiT = iT+off_Ash(iSym)
                  iiU = iU+off_Ash(jSym)
                  if (iiU > iiT) then
                    iiT = iU+off_Ash(jSym)
                    iiU = iT+off_Ash(iSym)
                  end if

                  iTU = iiU+nTri_Elem(iiT-1)
                  iiV = iV+off_Ash(kSym)
                  iiX = iX+off_Ash(lSym)
                  if (iiX > iiV) then
                    iiV = iX+off_Ash(lSym)
                    iiX = iV+off_Ash(kSym)
                  end if
                  iVX = iiX+nTri_Elem(iiV-1)
                  if (iVX > iTU) then
                    iTemp = iTU
                    iTU = iVX
                    iVX = iTemp
                  end if
                  iTUVX = iVX+nTri_Elem(iTU-1)
                  TUVX(iTUVX) = PUVX(iPUVX)
                end if

              end do
            end do
          end do
        end do
      end if

    end do
  end do
end do

return

end subroutine Get_TUVX
