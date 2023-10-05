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

subroutine CHO_rassi_twxy(irc,Scr,ChoV,TUVX,nAorb,JSYM,NUMV,DoReord)

use Cholesky, only: nSym
use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Data_Structures, only: SBA_Type, twxy_type
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: irc
integer(kind=iwp), intent(in) :: nAorb(*), JSYM, NUMV
type(twxy_type), intent(inout) :: Scr
type(SBA_Type), intent(in) :: ChoV
real(kind=wp), intent(_OUT_) :: TUVX(*)
logical(kind=iwp), intent(in) :: DoReord
integer(kind=iwp) :: iAorb(8), iRes, iSym, iSymt, iSymw, iSymx, iSymy, it, itG, itw, itwG, iw, iwG, ix, ixG, ixy, ixyG, iy, iyG, &
                     nTA, Ntw, Nxy

! ************************************************
if (NumV < 1) return

! Computing the integrals (TT|TT),(TW,TW) and (TW|XY)
!------------------------------------------------------
! (tw|xy)  <-  (tw|xy)  +  sum_J  L(tw,#J) * L(xy,#J)
!======================================================

do iSymy=1,nSym

  iSymx = Mul(iSymy,JSYM)

  Nxy = nAorb(iSymx)*nAorb(iSymy)

  if (Nxy > 0) then

    do iSymw=iSymy,nSym ! iSymw >= iSymy (particle symmetry)

      iSymt = Mul(iSymw,JSYM)

      Ntw = nAorb(iSymt)*nAorb(iSymw)

      if (Ntw > 0) then

        call DGEMM_('N','T',Ntw,Nxy,NumV,ONE,ChoV%SB(iSymt)%A3,Ntw,ChoV%SB(iSymx)%A3,Nxy,One,Scr%SB(iSymw,iSymy)%A,Ntw)

      end if

    end do

  end if

end do

! Reorder to the storage required by the RASSI program
!
! There is no permutational symmetry but only particle
! symmetry in the (tw|xy) integrals
!---------------------------------------------------------
if (DoReord) then

  iAorb(1) = 0
  do iSym=2,nSym
    iAorb(iSym) = iAorb(iSym-1)+nAorb(iSym-1)
  end do

  nTA = iAorb(nSym)+nAorb(nSym) ! total # active orbitals

  do iSymy=1,nSym

    iSymx = Mul(iSymy,JSYM)

    Nxy = nAorb(iSymx)*nAorb(iSymy)

    if (Nxy <= 0) cycle

    do iSymw=iSymy,nSym

      iSymt = Mul(iSymw,JSYM)

      Ntw = nAorb(iSymt)*nAorb(iSymw)

      if (Ntw <= 0) cycle

      do iy=1,nAorb(iSymy)

        iyG = iAorb(iSymy)+iy ! global index

        do ix=1,nAorb(iSymx)

          ixG = iAorb(iSymx)+ix

          ixy = ix+nAorb(iSymx)*(iy-1)

          ixyG = nTA*(iyG-1)+ixG ! global index

          do iw=1,nAorb(iSymw)

            iwG = iAorb(iSymw)+iw

            do it=1,nAorb(iSymt)

              itG = iAorb(iSymt)+it

              itw = it+nAorb(iSymt)*(iw-1)

              itwG = nTA*(iwG-1)+itG ! global index

              iRes = iTri(itwG,ixyG)

              TUVX(iRes) = Scr%SB(iSymw,iSymy)%A(itw,ixy)

            end do

          end do

        end do

      end do

    end do

  end do

end if

irc = 0

return

end subroutine CHO_rassi_twxy
