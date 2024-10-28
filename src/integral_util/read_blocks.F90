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

subroutine Read_Blocks(iTable,nBlocks,nBas,nIrrep,Buf,nBuf,iSO2Shell,nSOs,Bin,nBin,nQuad,G_Toc,iSO2cI,CutInt)

use Index_Functions, only: iTri, nTri_Elem
use SOAO_Info, only: iOffSO
use PSO_Stuff, only: Case_MP2, LuGam, LuGamma
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBlocks, iTable(6,nBlocks), nIrrep, nBas(0:nIrrep-1), nBuf, nSOs, iSO2Shell(nSOs), nBin, nQuad
real(kind=wp), intent(out) :: Buf(nBuf), Bin(2,nBin,nQuad), G_Toc(nQuad)
integer(kind=iwp), intent(out) :: iSO2cI(2,nSOs)
real(kind=wp), intent(in) :: CutInt
integer(kind=iwp) :: iAB, iAB_E, iAB_S, iAdr, iBlock, iBlockAdr, iBuf, iCD, iDisk, iiAB, iIrrep_A, iIrrep_B, iIrrep_C, iIrrep_D, &
                     Index_A, Index_AB, Index_ABCD, Index_B, Index_C, Index_CD, Index_D, iQuad, iShell_A, iShell_Ab, iShell_ABCD, &
                     iShell_B, iShell_C, iShell_CD, iShell_D, iSize, iSO_A_A, iSO_A_R, iSO_B_A, iSO_B_R, iSO_C_A, iSO_C_R, &
                     iSO_D_A, iSO_D_R, iType, iWrite, kBin, nA, nAB, nAB_Dist, nB, nC, nCD, nD, nDim_A, nDim_AB, nDim_B, nDim_C, &
                     nDim_CD, nDim_D
real(kind=wp) :: ABCD, BackChain
logical(kind=iwp) :: Triangular

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate table SO to contiguous index

!write(u6,*) 'nQuad=',nQuad
call Mk_SO2cI(iSO2cI,iSO2Shell,nSOs)
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate bins

Bin(1,nBin,:) = Zero  ! number of elements in bin
Bin(2,nBin,:) = -One  ! back chaining address
iDisk = 0 ! initial disk address
iWrite = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over symmetry blocks of gammas

!write(u6,*) 'iSO2Shell',iSO2Shell
!write(u6,*) 'iSO2cI(1)',(iSO2cI(1,i),i=1,nSOs)
!write(u6,*) 'iSO2cI(2)',(iSO2cI(2,i),i=1,nSOs)
iBlockAdr = 1
do iBlock=1,nBlocks
  iType = iTable(1,iBlock)
  iIrrep_A = iTable(2,iBlock)
  iIrrep_B = iTable(3,iBlock)
  iIrrep_C = iTable(4,iBlock)
  iIrrep_D = iTable(5,iBlock)
  !IND = iTable(6,iBlock)
  !write(u6,*) 'Irreps:',iIrrep_A,iIrrep_B,iIrrep_C,iIrrep_D

  nA = nBas(iIrrep_A)
  nB = nBas(iIrrep_B)
  nC = nBas(iIrrep_C)
  nD = nBas(iIrrep_D)
  !write(u6,*) 'nA,nB,nC,nD=',nA,nB,nC,nD
  if ((iType == 1) .or. (iType == 2)) then
    nAB = nTri_Elem(nA)
    nCD = nTri_Elem(nC)
    Triangular = .true.
  else
    nAB = nA*nB
    nCD = nC*nD
    Triangular = .false.
  end if
  if (nAB*nCD == 0) cycle

  nAB_dist = min(nAB,nBuf/nCD)
  iSO_A_r = 1
  iSO_B_r = 1
  do iiAB=1,nAB,nAB_dist
    iAB_s = iiAB
    iAB_e = min(nAB,iiAB+nAB_dist-1)

    ! This call is to an ACES2 routine. This has to be replaced once this
    ! code is used by MOLCAS. RL 2007-10-18.
    !           Call GetLst(Buf,iAB_s,nAB_dist,2,iType,IND)
    !*******************************************************************
    ! Jonas B 2010. This code is used for reading nonseparable two-electron
    !               density matrices when doing conventional mp2-gradients.
    if (Case_mp2) then
      iSize = min(nAB_dist*nCD,(nAB-(iiAB-1))*nCD)
      iAdr = iBlockAdr+nCD*(iiAB-1)
      call dDaFile(LuGam,2,Buf,iSize,iAdr)
    end if

    !*******************************************************************
    !call RecPrt('Aces_Gamma: from GetLst',' ',Buf,iAB_e-iAB_s+1,nCD)

    iBuf = 0
    do iAB=iAB_s,iAB_e

      iSO_A_a = iOffSO(iIrrep_A)+iSO_A_r
      iSO_B_a = iOffSO(iIrrep_B)+iSO_B_r
      iShell_A = iSO2Shell(iSO_A_a)
      iShell_B = iSO2Shell(iSO_B_a)
      iShell_AB = iTri(iShell_A,iShell_B)

      Index_A = iSO2cI(1,iSO_A_a)
      Index_B = iSO2cI(1,iSO_B_a)
      nDim_A = iSO2cI(2,iSO_A_a)
      nDim_B = iSO2cI(2,iSO_B_a)
      nDim_AB = nDim_A*nDim_B
      if (iShell_A > iShell_B) then
        Index_AB = (Index_B-1)*nDim_A+Index_A
      else if (iShell_A == iShell_B) then
        Index_AB = iTri(Index_A,Index_B)
      else
        Index_AB = (Index_A-1)*nDim_B+Index_B
      end if

      iSO_C_r = 1
      iSO_D_r = 1
      do iCD=1,nCD

        iBuf = iBuf+1
        !iBuf = (iCD-1)*nAB_dist+iAB
        ABCD = Buf(iBuf)

        ! Jonas B 2010. If calculating mp2-integrals the cutoff for small
        ! gammas should not be used at the moment.
        if ((.not. Case_mp2) .and. (abs(ABCD) < CutInt)) cycle
        iSO_C_a = iOffSO(iIrrep_C)+iSO_C_r
        iSO_D_a = iOffSO(iIrrep_D)+iSO_D_r
        iShell_C = iSO2Shell(iSO_C_a)
        iShell_D = iSO2Shell(iSO_D_a)
        iShell_CD = iTri(iShell_C,iShell_D)
        iShell_ABCD = iTri(iShell_AB,iShell_CD)

        ! Compute canonical compound index

        Index_C = iSO2cI(1,iSO_C_a)
        Index_D = iSO2cI(1,iSO_D_a)
        nDim_C = iSO2cI(2,iSO_C_a)
        nDim_D = iSO2cI(2,iSO_D_a)
        nDim_CD = nDim_C*nDim_D
        if (iShell_C > iShell_D) then
          Index_CD = (Index_D-1)*nDim_C+Index_C
        else if (iShell_C == iShell_D) then
          Index_CD = iTri(Index_C,Index_D)
        else
          Index_CD = (Index_C-1)*nDim_D+Index_D
        end if
        if (iShell_AB > iShell_CD) then
          Index_ABCD = (Index_CD-1)*nDim_AB+Index_AB
        else if (iShell_AB == iShell_CD) then
          Index_ABCD = iTri(Index_AB,Index_CD)
        else
          Index_ABCD = (Index_AB-1)*nDim_CD+Index_CD
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Store information in the appropriate bin

        !write(u6,*) 'ABCD, Index_ABCD=',ABCD, Index_ABCD
        !write(u6,*) 'iSO:',iSO_A_a,iSO_B_a,iSO_C_a,iSO_D_a
        !write(u6,*) 'iShell:',iShell_A,iShell_B,iShell_C,iShell_D
        !write(u6,*) 'cind:',Index_A,Index_B,Index_C,Index_D
        !write(u6,*) 'ndims:',nDim_A,nDim_B,nDim_C,nDim_D
        !write(u6,*) 'iShell_AB,iShell_CD=',iShell_AB,iShell_CD
        !write(u6,*) 'nDim_AB,nDim_CD=',nDim_AB,nDim_CD
        !write(u6,*) 'iShell_ABCD=',iShell_ABCD
        kBin = int(Bin(1,nBin,iShell_ABCD))+1
        Bin(1,kBin,iShell_ABCD) = ABCD
        Bin(2,kBin,iShell_ABCD) = real(Index_ABCD,kind=wp)
        Bin(1,nBin,iShell_ABCD) = real(kBin,kind=wp)  ! Update counter
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Write bin to disk if full.

        if (kBin == nBin-1) then
          BackChain = real(iDisk,kind=wp)
          call dDaFile(LuGamma,iWrite,Bin(:,:,iShell_ABCD),2*nBin,iDisk)
          Bin(1,nBin,iShell_ABCD) = Zero
          Bin(2,nBin,iShell_ABCD) = BackChain
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Increment indices for C and D

        if (.not. Triangular) then
          !iSO_C_r = iSO_C_r+1
          !if (iSO_C_r > nBas(iIrrep_C)) then
          !  iSO_C_r = 1
          !  iSO_D_r = iSO_D_r+1
          !end if
          iSO_D_r = iSO_D_r+1
          if (iSO_D_r > nBas(iIrrep_D)) then
            iSO_D_r = 1
            iSO_C_r = iSO_C_r+1
          end if
        else
          ! Lower triangular
          iSO_D_r = iSO_D_r+1
          if (iSO_D_r > iSO_C_r) then
            iSO_C_r = iSO_C_r+1
            iSO_D_r = 1
          end if
          ! Upper triangular
          !iSO_C_r = iSO_C_r+1
          !if (iSO_C_r > nBas(iIrrep_C)) then
          !  iSO_D_r = iSO_D_r+1
          !  iSO_C_r = iSO_D_r
          !end if
        end if

      end do  ! iCD
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Increment indices for A and B

      if (.not. Triangular) then
        !iSO_A_r = iSO_A_r+1
        !if (iSO_A_r > nBas(iIrrep_A)) then
        !  iSO_A_r = 1
        !  iSO_B_r = iSO_B_r+1
        !end if
        iSO_B_r = iSO_B_r+1
        if (iSO_B_r > nBas(iIrrep_B)) then
          iSO_B_r = 1
          iSO_A_r = iSO_A_r+1
        end if
      else
        ! Lower triangular
        iSO_B_r = iSO_B_r+1
        if (iSO_B_r > iSO_A_r) then
          iSO_A_r = iSO_A_r+1
          iSO_B_r = 1
        end if
        ! Upper triangular
        !iSO_A_r = iSO_A_r+1
        !if (iSO_A_r > nBas(iIrrep_A)) then
        !  iSO_B_r = iSO_B_r+1
        !  iSO_A_r = iSO_B_r
        !end if
      end if

    end do     ! iAB

  end do       ! iiAB
  iBlockAdr = iBlockAdr+nAB*nCD

end do         ! iBlock
!                                                                      *
!***********************************************************************
!                                                                      *
! Flush the bins to disk and save disk addresses

!write(u6,*) 'Flush bins!'
do iQuad=1,nQuad
  BackChain = real(iDisk,kind=wp)
  call dDaFile(LuGamma,iWrite,Bin(:,:,iQuad),2*nBin,iDisk)
  !lxx = Int(Bin(1,nBin,iQuad))
  !write(u6,*) 'lxx=',lxx
  !call RecPrt('Bins',' ',Bin(1,1,iQuad),2,lxx)
  G_Toc(iQuad) = BackChain
end do

! Observe that backchain list is only stored in core!
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Read_Blocks
