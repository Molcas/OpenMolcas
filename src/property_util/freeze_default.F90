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

subroutine Freeze_Default(iANr,nShell,lMax)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iANr, lMax
integer(kind=iwp), intent(out) :: nShell(0:lMax)
integer(kind=iwp) :: i
integer(kind=iwp), parameter :: nAtoms = 106, &
                                iDefaults(0:3,0:nAtoms) = reshape([ &
                                                                   ! First row, H-He (and ghost): no frozen orbitals
                                                                   0,0,0,0, &
                                                                   0,0,0,0, &
                                                                   0,0,0,0, &
                                                                   ! Second row, Li-Ne: 1s
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   1,0,0,0, &
                                                                   ! Third row, Na-Si: 1s2s
                                                                   2,0,0,0, &
                                                                   2,0,0,0, &
                                                                   2,0,0,0, &
                                                                   2,0,0,0, &
                                                                   ! Third row, P-Ar: 1s2s2p
                                                                   2,1,0,0, &
                                                                   2,1,0,0, &
                                                                   2,1,0,0, &
                                                                   2,1,0,0, &
                                                                   ! Fourth row, K-Cu: 1s2s2p3s
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   3,1,0,0, &
                                                                   ! Fourth row, Zn-Kr: 1s2s2p3s3p3d
                                                                   3,2,1,0, &
                                                                   3,2,1,0, &
                                                                   3,2,1,0, &
                                                                   3,2,1,0, &
                                                                   3,2,1,0, &
                                                                   3,2,1,0, &
                                                                   ! Fifth row, Rb-Ag: 1s2s2p3s3p4s3d
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   4,2,1,0, &
                                                                   ! Fifth row, Cd-Xe: 1s2s2p3s3p4s3d4p
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   4,3,1,0, &
                                                                   ! Sixth row, Cs-Yb: 1s2s2p3s3p4s3d4p5s4d5p
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   4,3,2,0, &
                                                                   ! Sixth row, Lu-Au: 1s2s2p3s3p4s3d4p4d4f
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   4,3,2,1, &
                                                                   ! Sixth row, Hg-Rn: 1s2s2p3s3p4s3d4p4d4f5s5p
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   5,4,2,1, &
                                                                   ! Seventh row, Fr-Cm: 1s2s2p3s3p4s3d4p5s4d5p4f5d
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   5,4,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,1, &
                                                                   6,5,3,2, &
                                                                   6,5,3,2, &
                                                                   6,5,3,2  &
                                                                  ],shape(iDefaults))

if ((iANr < 0) .or. (iANr > nAtoms)) then
  write(u6,*) 'Freeze_Defaults: iAnr is out of range!'
  write(u6,*) 'iANr=',iANr
  call Abend()
end if

nShell(:) = 0
do i=0,min(lMax,3)
  nShell(i) = iDefaults(i,iAnr)
end do

return

end subroutine Freeze_Default
