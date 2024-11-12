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
! Copyright (C) 2004, Par Soderhjelm                                   *
!***********************************************************************

subroutine XFMoment(lMax,Cavxyz,Tmom,nCavxyz_,Org)
!***********************************************************************
!                                                                      *
!     Object:  Calculate total moment of XF multipoles, add to Cavxyz  *
!                                                                      *
!     Authors: P. Soderhjelm                                           *
!              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
!                                                                      *
!              November 2004                                           *
!***********************************************************************

use Index_Functions, only: nTri3_Elem1
use External_Centers, only: nOrd_XF, nXF, XF
use Phase_Info, only: iPhase
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lmax, nCavxyz_
real(kind=wp), intent(inout) :: Cavxyz(nCavxyz_)
real(kind=wp), intent(out) :: Tmom(nCavxyz_), Org(3)
integer(kind=iwp) :: i, iChAtm, iChxyz, iDum, iStb(0:7), j, jCoSet(0:7,0:7), nInp, nStb
real(kind=wp) :: A(3), Tco(3)

if (nOrd_XF < 0) return
if (nOrd_XF > lMax) then
  call WarningMessage(2,'nOrd_XF > lMax')
  call Abend()
end if
nInp = nTri3_Elem1(nOrd_XF)
Org(:) = Zero
do i=1,nXF

  ! Generate Stabilizer of C

  ! IFG: "A" was undefined, is this the right point?
  A(1:3) = XF(1:3,i)
  iChxyz = iChAtm(A)
  iDum = 0
  call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)
  do j=0,nIrrep/nStb-1
    Tmom(:) = Zero
    Tmom(1:nInp) = XF(4:3+nInp,i)
    TCo(1:3) = XF(1:3,i)
    Tco(1) = Tco(1)*real(iPhase(1,jCoSet(j,0)),kind=wp)
    Tco(2) = Tco(2)*real(iPhase(2,jCoSet(j,0)),kind=wp)
    Tco(3) = Tco(3)*real(iPhase(3,jCoSet(j,0)),kind=wp)
    if (nOrd_XF > 0) then
      Tmom(2) = Tmom(2)*real(iPhase(1,jCoSet(j,0)),kind=wp)   !Dx
      Tmom(3) = Tmom(3)*real(iPhase(2,jCoSet(j,0)),kind=wp)   !Dy
      Tmom(4) = Tmom(4)*real(iPhase(3,jCoSet(j,0)),kind=wp)   !Dz
      if (nOrd_XF > 1) then
        Tmom(6) = Tmom(6)*real(iPhase(1,jCoSet(j,0))*iPhase(2,jCoSet(j,0)),kind=wp) !Qxy
        Tmom(7) = Tmom(7)*real(iPhase(1,jCoSet(j,0))*iPhase(3,jCoSet(j,0)),kind=wp) !Qxz
        Tmom(9) = Tmom(9)*real(iPhase(2,jCoSet(j,0))*iPhase(3,jCoSet(j,0)),kind=wp) !Qyz
      end if
    end if
    call ReExpand(Tmom,1,nCavxyz_,Tco,Org,1,lMax)
    Cavxyz(:) = Cavxyz(:)+Tmom(:)
  end do
end do

end subroutine XFMoment
