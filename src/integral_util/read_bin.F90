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

subroutine Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Gmma,nGamma,LuGamma,Bin,lBin)

use Index_Functions, only: iTri
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iShell_A, iShell_B, iShell_C, iShell_D, nQuad, nGamma, LuGamma, lBin
real(kind=wp), intent(in) :: G_Toc(nQuad)
real(kind=wp), intent(out) :: Gmma(nGamma), Bin(2,lBin)
integer(kind=iwp) :: iDisk, iGamma, iShell_AB, iShell_ABCD, iShell_CD, jGamma, lGamma
integer(kind=iwp), parameter :: iRead = 2

Gmma(:) = Zero

iShell_AB = iTri(iShell_A,iShell_B)
iShell_CD = iTri(iShell_C,iShell_D)
iShell_ABCD = iTri(iShell_AB,iShell_CD)

!write(u6,*) 'Reading Gammas for shell quadruplet ',iShell_ABCD

iDisk = int(G_Toc(iShell_ABCD))
do while (iDisk >= 0)
  call dDaFile(LuGamma,iRead,Bin,2*lBin,iDisk)
  lGamma = int(Bin(1,lBin))
  iDisk = int(Bin(2,lBin))

  do iGamma=1,lGamma
    jGamma = int(Bin(2,iGamma))
    if (jGamma > nGamma) then
      call WarningMessage(2,'Read_Bin: jGamma > nGamma')
      call Abend()
    end if
    Gmma(jGamma) = Bin(1,iGamma)
    !write(u6,*) Gmma(jGamma),jGamma
  end do

end do

end subroutine Read_Bin
