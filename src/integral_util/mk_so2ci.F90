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

subroutine Mk_SO2cI(SO2cI,iSO2Shell,nSOs)

implicit none
integer nSOs
integer SO2cI(2,nSOs), iSO2Shell(nSOs)
integer nShell, iShell, Index, iSO

!                                                                      *
!***********************************************************************
!                                                                      *
! Set up table SO to contigues index over the shell

!write(6,*) 'Enter SO2cI'
call ICopy(2*nSOs,[0],0,SO2cI,1)
call Nr_Shells(nShell)
do iShell=1,nShell

  ! Generate contigues index for this shell

  Index = 0
  do iSO=1,nSOs
    if (iSO2Shell(iSO) == iShell) then
      Index = Index+1
      SO2cI(1,iSO) = Index
    end if
  end do

  ! Store dimension for this shell

  do iSO=1,nSOs
    if (iSO2Shell(iSO) == iShell) SO2cI(2,iSO) = Index
  end do

end do ! iShell
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_SO2cI
