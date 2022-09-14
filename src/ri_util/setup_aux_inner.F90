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

subroutine Setup_Aux_Inner(iSOShl,nSO,iShlSO,nBasSh,nShell,nIrrep,nBas,iSSOff,nij_Shell,iShij,nBas_Aux,nChV,iTOffs)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSO, iSOShl(nSO), nShell, nIrrep, nBas(0:nIrrep-1), nij_Shell, iShij(2,nij_Shell), &
                                 nBas_Aux(0:nIrrep-1), nChV(0:nIrrep-1)
integer(kind=iwp), intent(out) :: iShlSO(nSO), nBasSh(0:nIrrep-1,nShell), iSSOff(0:nIrrep-1,0:nIrrep-1,nij_Shell), &
                                  iTOffs(3,0:nIrrep-1)
integer(kind=iwp) :: iAcc, iBas, iIrrep, ijIrrep, ijShell, iOff_V12, iShell, iSO, iSO_Shl, iTtmp(0:7), jIrrep, jShell, nA, nab, &
                     nAux, nB, nI

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate index array for relative index within the shell and irrep

!call iVcPrt('iSOShl',' ',iSOShl,nSO)
iSO = 0
do iIrrep=0,nIrrep-1
  do iShell=1,nShell

    iSO_Shl = 0
    do iBas=iSO+1,iSO+nBas(iIrrep)
      if (iSOShl(iBas) == iShell) then
        iSO_Shl = iSO_Shl+1
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Save the relative index within the shell and irrep
        ! of a given absolute SO index.

        iShlSO(iBas) = iSO_Shl
      end if
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Save the total number of basis functions a specific shell
    ! has in a given irrep.

    nBasSh(iIrrep,iShell) = iSO_Shl
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do
  iSO = iSO+nBas(iIrrep)
end do
!                                                                      *
!***********************************************************************
!
! Initialize

iTOffs(:,:) = 0
iSSOff(:,:,:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute offsets within the symmetry block for a fixed shell pair.

! Note that for each pair of valence shells all the products which
! are of the same irrep are consecutive.

do ijShell=1,nij_Shell
  iShell = iShij(1,ijShell)
  jShell = iShij(2,ijShell)
  iTtmp(0:nIrrep-1) = 0
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iShell > jShell) then   ! iShell > jShell

    do jIrrep=0,nIrrep-1
      nB = nBasSh(jIrrep,jShell)
      do iIrrep=0,nIrrep-1
        nA = nBasSh(iIrrep,iShell)

        ijIrrep = Mul(iIrrep+1,jIrrep+1)-1
        iSSOff(iIrrep,jIrrep,ijShell) = iTtmp(ijIrrep)
        nab = na*nb
        iTtmp(ijIrrep) = iTtmp(ijIrrep)+nab

      end do
    end do

  else       ! iShell = jShell

    do iIrrep=0,nIrrep-1
      nA = nBasSh(iIrrep,iShell)
      do jIrrep=0,iIrrep
        nB = nBasSh(jIrrep,jShell)

        ijIrrep = Mul(iIrrep+1,jIrrep+1)-1
        iSSOff(iIrrep,jIrrep,ijShell) = iTtmp(ijIrrep)
        iSSOff(jIrrep,iIrrep,ijShell) = iTtmp(ijIrrep)
        nab = na*nb
        if (iIrrep == jIrrep) nab = nTri_Elem(na)
        iTtmp(ijIrrep) = iTtmp(ijIrrep)+nab
      end do
    end do

  end if

  do iIrrep=0,nIrrep-1
    iTOffs(3,iIrrep) = iTOffs(3,iIrrep)+iTtmp(iIrrep)
  end do

  ! Now update the index to be the total offset within a slice
  ! for a fixed shell-pair

  iAcc = 0
  do ijIrrep=0,nIrrep-1
    do iIrrep=0,nIrrep-1
      jIrrep = Mul(ijIrrep+1,iIrrep+1)-1
      iSSOff(iIrrep,jIrrep,ijShell) = iSSOff(iIrrep,jIrrep,ijShell)+iAcc
    end do
    nI = nBas_Aux(ijIrrep)
    if (ijIrrep == 0) nI = nI-1
    iAcc = iAcc+nI*iTtmp(ijIrrep)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'iSSOff'
write(u6,*)
do ijShell=1,nij_Shell
  iShell = iShij(1,ijShell)
  jShell = iShij(2,ijShell)
  write(u6,*)
  write(u6,*) 'iShell,jShell=',iShell,jShell
  write(u6,*)
  do i=0,nIrrep-1
    write(u6,'(8I4)') (iSSOff(i,j,ijShell),j=0,nIrrep-1)
  end do
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up pointers for the J12 matrix and compute total size of the
! 3-center integrals.

iOff_V12 = 0
do iIrrep=0,nIrrep-1
  iTOffs(1,iIrrep) = nChV(iIrrep) ! # of vectors
  nAux = nBas_Aux(iIrrep)
  if (iIrrep == 0) nAux = nAux-1
  iTOffs(2,iIrrep) = iOff_V12
  iOff_V12 = iOff_V12+nAux**2
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' iSO, iShlSO(iSO), relative index in irrep'
do jSO=1,iSO
  write(u6,*) jSO,iShlSO(jSO)
end do

write(u6,*)
write(u6,*) ' iShell: number of basis functions in each irrep'
do iShell=1,nShell
  write(u6,*) iShell,':',(nBasSh(iIrrep,iShell),iIrrep=0,nIrrep-1)
end do
write(u6,*)
write(u6,*) 'iSSOff'
write(u6,*)
do ijShell=1,nij_Shell
  iShell = iShij(1,ijShell)
  jShell = iShij(2,ijShell)
  write(u6,*)
  write(u6,*) 'iShell,jShell=',iShell,jShell
  write(u6,*)
  do i=0,nIrrep-1
    write(u6,'(8I4)') (iSSOff(i,j,ijShell),j=0,nIrrep-1)
  end do
end do
write(u6,*)
write(u6,*) 'iTOffs'
write(u6,*)
do i=0,nIrrep-1
  write(u6,*) (iTOffs(j,i),j=1,3)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Setup_Aux_Inner
