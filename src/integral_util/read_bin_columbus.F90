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

subroutine Read_Bin_Columbus(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Gamma,nGamma,LuGamma,Bin,lBin)

use Constants, only: Zero

implicit none
#include "SysDef.fh"
integer iShell_A, iShell_B, iShell_C, iShell_D, nQuad, nGamma, LuGamma, lBin
real*8 G_Toc(nQuad), Bin(2,lBin), gamma(nGamma)
integer iShell_AB, iShell_CD, iShell_ABCD
integer, parameter :: iRead = 2
integer idisk_save
integer iDisk, lGamma, iGamma, jGamma
integer i, j, iTri
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

gamma(:) = Zero

iShell_AB = iTri(iShell_A,iShell_B)
iShell_CD = iTri(iShell_C,iShell_D)
iShell_ABCD = iTri(iShell_AB,iShell_CD)

!write(u6,*) 'Reading Gammas for shell quadruplet ',iShell_ABCD

iDisk = int(G_Toc(iShell_ABCD))
do while (iDisk >= 0)
  idisk_save = idisk
  call dDaFile(LuGamma,iRead,Bin,2,iDisk)
  lGamma = int(Bin(1,1))
  if (lGamma > lBin) then
    call WarningMessage(2,'Read_Bin_Columbus: lGamma > lbin')
    call Abend()
  end if
  iDisk = idisk_save
  call dDaFile(LuGamma,iRead,Bin,2*lGamma,iDisk)
  iDisk = int(Bin(2,1))

  do iGamma=2,lGamma
    jGamma = int(Bin(2,iGamma))
    if (jGamma > nGamma) then
      call WarningMessage(2,'Read_Bin_Columbus: jGamma > nGamma')
      call Abend()
    end if
    gamma(jGamma) = Bin(1,iGamma)
    !write(u6,*) Gamma(jGamma),jGamma
  end do

end do

end subroutine Read_Bin_Columbus
