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
! Copyright (C) 2009, Giovanni Li Manni                                *
!               2009, Francesco Aquilante                              *
!***********************************************************************

subroutine splitCTL(LW1,TUVX,IFINAL,iErrSplit)
!***********************************************************************
!     CI Hamiltonian Matrix elements reader                            *
!     calling arguments:                                               *
!     LW1     : Memory pointer to active Fock matrix                   *
!               array of real                                          *
!     TUVX    : array of real                                          *
!               two-electron integrals (tu!vx)                         *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Written by:                                                      *
!     G. Li Manni (GLMJ) and F. Aquilante                              *
!     University of Geneva, Switzerland, 2009                          *
!                                                                      *
!***********************************************************************

use csfbas, only: CONF, NAEL, NBEL
use timers, only: C_Dress, C_get_Cm, W_Dress, W_get_Cm
use lucia_data, only: DFTP, DTOC, CFTP, IREOTS
use splitcas_data, only: EnInSplit, EnerSplit, FordSplit, gapSpli, iDimBlockA, iDimBlockACNF, iterSplit, lRootSplit, MxIterSplit, &
                         NumSplit, percSpli, PerSplit, ThrSplit
use rasscf_global, only: EMY, ENER, ExFac, IADR15, ICICH, iCIOnly, ICIRST, ITER, ITERCI, iTOC, n_Keep, NAC
use general_data, only: JOBIPH, JOBOLD, LUDAVID, NACTEL, NCONF, STSYM
use spinfo, only: NCNASM, NDET
use output_ras, only: IPRLOC
use printlevel, only: DEBUG
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, auToeV
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: LW1(*), TUVX(*)
integer(kind=iwp), intent(in) :: IFINAL
integer(kind=iwp), intent(out) :: iErrSplit
integer(kind=iwp) :: iCaseSplit, iDimBlockTri, iDisk, idx, iJOB, iPrint, IPRLEV, j, k, MXSpli, MXXWS, nAAblock
real(kind=wp) :: C_ABlockDim_sel(2), condition, CSplitTot(2), Ctime(2), diffSplit, ECORE, EnFinSplit, SpliNor, W_ABlockDim_sel(2), &
                 WSplitTot(2), Wtime(2)
character(len=80) :: String
logical(kind=iwp) :: DBG, Exists
integer(kind=iwp), allocatable :: IPCNF(:), IPCNFtot(:), IPCSFtot(:), iSel(:), vkcnf(:)
real(kind=wp), allocatable :: AABlock(:), CIVEC(:), DHAM(:), Diag(:), DiagCNF(:), HONE(:,:), Scr(:), SplitE(:), SplitV(:,:), &
                              Tmp1(:), Tmp2(:), TotSplitV(:)
real(kind=wp), external :: ddot_

#include "macros.fh"
unused_var(IFINAL)

IPRLEV = IPRLOC(3)
IPRINT = IPRLEV

DBG = IPRLEV >= DEBUG

iErrSplit = 0
!-----------------------------------------------------------------------
!    INITIALIZE THE DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------
!
!1) In posizione 1 nella seguente routine ho messo 'lRootSplit' invece
!   di 'lRoots', in quanto 'lRoots' viene specificato nell'input del
!   calcolo con la keyword CIROOT (secondo numero intero scritto nella
!   linea sotto CIROOT) ma sottomettendo un calcolo SplitCAS non andiamo
!   ad usare tale keyword (almeno non per il momento!).
!
!2) In posizione 4 nella seguente routine ho messo 'nConf' invece di 'nSel'
!   rispetto alla routine originale in 'davctl'.
!
! Qui nConf rappresenta il numero totale di CSFs.

call cwtime(CSplitTot(1),WSplitTot(1))
!call Ini_David(lRootSplit,nConf,nDet,nconf,nAc,LuDavid)
call Ini_David(1,nConf,nDet,nconf,n_keep,nAc,LuDavid)
!-----------------------------------------------------------------------
!     COMPUTE THE DIAGONAL ELEMENTS OF THE COMPLETE HAMILTONIAN
!-----------------------------------------------------------------------
! CIVEC: TEMPORARY CI VECTOR IN CSF BASIS
call mma_allocate(CIVEC,NCONF,label='CIVEC')

!if (IFINAL == 2) then ! to avoid last diagonalization
!  CIVEC(:) = Zero
!  !call Load_tmp_CI_vec(1,1,nConf,CIVEC,LuDavid)
!  call Load_CI_vec(1,nConf,CIVEC,LuDavid)
!  !call dDaFile(JOBIPH,2,CIVEC,nConf,LuDavid)
!  if (DBG) then
!    write(u6,*) 'LuDavid',LuDavid
!    write(String,'(A)') 'Final=2 : CI-coeff in SplitCAS'
!    call dVcPrt(String,' ',CIVEC,nConf)
!  end if
!
!  if (NAC == 0) then
!    ENER(1,ITER) = EMY
!  else
!    !ENER(lRootSplit,ITER) = ENER(lRootSplit,ITER-1)
!    if (DBG) then
!      write(u6,*) 'lRootSplit :',lRootSplit
!      write(u6,*) 'ITER :',ITER
!      write(u6,*) 'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
!    end if
!  end if
!  call mma_deallocate(CIVEC)
!  return
!end if

if (NAC > 0) call CIDIA(NCONF,STSYM,CIVEC,LUDAVID)
!***********************************************************************
! iCaseSplit = 1  : there is NOT CI-RESTART.
! iCaseSplit = 2  : there is CIRESTART. The code will read the CI
!                   coefficients from the JOBOLD file.
!***********************************************************************
iCaseSplit = 1
if (ICIRST /= 0) iCaseSplit = 2

if (iCaseSplit == 1) then ! There is NO CIRST
  if (iDimBlockA /= nconf) then ! AA Block is smaller than the full Hamiltonian matrix, then SplitCAS will be performed.

    call mma_allocate(HONE,NAC,NAC,label='HONE')
    ! EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE
    call SQUARE(LW1,HONE,NAC,1,NAC)
    !call mma_allocate(IPCNF,NCNASM(STSYM),label='IPCNF')

    call cwtime(C_ABlockDim_sel(1),W_ABlockDim_sel(1))
    ECORE = Zero
    if (NumSplit) then
      !*****************************************************************
      ! In such case the splitting is done according to a NUMERICAL    *
      ! selection.                                                     *
      ! In ipCSFSplit the following quantities are calculated:         *
      !    1) The index array for the CSFs and CNFs                    *
      !    2) The final values of iDimBlockA and iDimBlockACNF         *
      !    3) The diagonal array of the Hamiltonian matrix             *
      !*****************************************************************
      call mma_allocate(Diag,nConf,label='DiagCSF')
      call mma_allocate(IPCSFtot,nConf,label='IPCSFtot')
      call mma_allocate(IPCNFtot,NCNASM(STSYM),label='IPCNFtot')
      call mma_maxDBLE(MXXWS)
      call mma_allocate(Scr,MXXWS,label='EXHSCR')
      MXSpli = iDimBlockA
      !nAAblock = MXSpli*(MXSpli+1)/2
      call ipCSFSplit(Diag,IPCSFtot,IPCNFtot,nConf,MXSpli,DTOC,DFTP,CONF,STSYM,HONE,ECORE,NAC,Scr,NCNASM(STSYM), &
                      NAEL+NBEL,NAEL,NBEL,CIVEC,TUVX,IPRINT,ExFac,IREOTS)
      !call DVCPRT('Diagonal elements of Hamilt. matrix in CSF',' ',Diag,nConf)
      call mma_deallocate(Scr)
      !call mma_deallocate(Diag)
      !nAAblock = iDimBlockA*(iDimBlockA+1)/2
    end if

    if (EnerSplit .or. PerSplit) then
      !*****************************************************************
      ! COMPUTE the index array for THE DIAGONAL ELEMENTS OF THE       *
      ! COMPLETE HAMILTONIAN in energetic order.                       *
      !*****************************************************************
      call mma_allocate(Diag,nConf,label='DiagCSF')
      call mma_allocate(DiagCNF,NCNASM(STSYM),label='DiagCNF')
      call mma_allocate(IPCSFtot,nConf,label='IPCSFtot')
      call mma_allocate(IPCNFtot,NCNASM(STSYM),label='IPCNFtot')
      call mma_maxDBLE(MXXWS)
      call mma_allocate(Scr,MXXWS,label='EXHSCR')
      ! 'GapSpli' comes from the input in eV
      ! 'condition' goes to DiagOrd in hartree if EnerSplit
      ! 'condition' goes to DiagOrd as a percentage if PerSplit
      if (EnerSplit) condition = gapSpli/auToeV
      if (PerSplit) condition = percSpli
      call DiagOrd(Diag,DiagCNF,IPCSFtot,IPCNFtot,nConf,condition,ITER,DTOC,DFTP,CONF,STSYM,HONE,ECORE,NAC,Scr, &
                   NCNASM(STSYM),(NAEL+NBEL),NAEL,NBEL,TUVX,IPRINT,ExFac,IREOTS)
      if (DBG) then
        call DVCPRT('Diagonal elements of Hamilt. matrix in CSF',' ',Diag,nConf)
        call DVCPRT('Diagonal elements of Hamilt. matrix in CNF',' ',DiagCNF,NCNASM(STSYM))
        call IVCPRT('Index Array in CSF',' ',IPCSFtot,nConf)
        call IVCPRT('Index Array in CNF',' ',IPCNFtot,NCNASM(STSYM))
        write(u6,*) 'iDimBlockACNF from DiagOrd : ',iDimBlockACNF
        write(u6,*) 'iDimBlockA from DiagOrd : ',iDimBlockA
        call xflush(u6)
      end if
      call mma_deallocate(Scr)
      !call mma_deallocate(Diag)
      call mma_deallocate(DiagCNF)
      !call mma_deallocate(IPCNFtot)
    end if

    if (DBG) then
      call cwtime(C_ABlockDim_sel(2),W_ABlockDim_sel(2))
      write(u6,*) 'Time needed to select CSFs in A-Block:'
      write(u6,*) 'CPU timing : ',C_ABlockDim_sel(2)-C_ABlockDim_sel(1)
      write(u6,*) 'W. timing  : ',W_ABlockDim_sel(2)-W_ABlockDim_sel(1)
    end if
    !*******************************************************************
    ! Let's start with SPLITCAS code                                   *
    !*******************************************************************
    !call cwtime(C_iterative1,W_iterative1)
    call chksplit()
    !*******************************************************************
    if (DBG) then
      write(u6,*) 'iDimBlockA after selection : ',iDimBlockA
      write(u6,*) 'iDimBlockACNF after selection : ',iDimBlockACNF
    end if
    ! calculate the dressed HAMILTONIAN matrix (STEP 2)
    iDimBlockTri = iDimBlockA*(iDimBlockA+1)/2
    !iBCSF = nConf-iDimBlockA
    !iBCNF = NCNASM(STSYM)-NPCNF
    !iBMAX = MAX(iBCSF,iDimBlockA)
    call mma_allocate(AABlock,iDimBlockTri,label='AAblock')
    call mma_allocate(DHAM,iDimBlockTri,label='DHAM')
    !call mma_allocate(BVEC,iBMAX,label='BVEC')
    AABlock(:) = Zero
    if (iter == 1) then
      EnInSplit = Diag(lrootSplit)
    end if
    if (DBG) then
      write(u6,*) 'Initial Energy in SplitCAS : ',EnInSplit
    end if
    diffSplit = ThrSplit+One

    iterSplit = 0
    call cwtime(Ctime(1),Wtime(1))
    call mma_allocate(SplitE,iDimBlockA,label='SplitE')
    call mma_allocate(SplitV,iDimBlockA,iDimBlockA,label='SplitV')
    do while ((diffSplit > ThrSplit) .and. (iterSplit < MxIterSplit))
      iterSplit = iterSplit+1
      if (DBG) then
        write(u6,*) '*************** Iteration SplitCAS =',iterSplit
      end if
      !call Compute_Umn(BVEC,NPCNF,NCNASM(STSYM),EnInSplit,NPCNF+1,1,DHAM)
      !call SPLITCSF(AABlock,EnInSplit,DHAM,
      call get_Umn(AABlock,EnInSplit,DHAM,IPCSFtot,IPCNFtot,nconf,DTOC,DFTP,CONF,STSYM,HONE,ECORE,NAC, &
                   NCNASM(STSYM),NAEL+NBEL,NAEL,NBEL,iDimBlockA,iDimBlockACNF,TUVX,iterSplit,ITER,IPRINT,ExFac,IREOTS)
      call xflush(u6)
      if (DBG) then
        call TRIPRT('AA block of the Hamiltonian Matrix',' ',AABlock,iDimBlockA)
        call TRIPRT('Dressed AA block Hamiltonian Matrix',' ',DHAM,iDimBlockA)
        call xflush(u6)
      end if
      !call mma_deallocate(BVEC)
      !*****************************************************************
      ! Dressed Hamiltonian DIAGONALIZATION                            *
      !*****************************************************************
      call unitmat(SplitV,iDimBlockA)
      call NIdiag(DHAM,SplitV,iDimBlockA,iDimBlockA)
      call JACORD(DHAM,SplitV,iDimBlockA,iDimBlockA)
      do idx=1,iDimBlockA
        SplitE(idx) = DHAM(idx*(idx+1)/2)
      end do
      EnFinSplit = SplitE(lrootSplit)
      if (DBG) then
        call IVCPRT('CSFs included ',' ',IPCSFtot,iDimBlockA)
        call RECPRT('Eigenvec. of dressed Hamiltonian',' ',SplitV,iDimBlockA,iDimBlockA)
        call DVCPRT('Eigenval. of dressed Hamiltonian',' ',SplitE,iDimBlockA)
      end if
      diffSplit = abs(EnFinSplit-EnInSplit)
      if (DBG) write(u6,*) 'Energy diff in SplitCAS :',diffSplit
      if (DBG) then
        if (iterSplit == 1) write(u6,*) 'IterSplit   Final Energy in SplitCAS'
        write(u6,*) IterSplit,EnFinSplit
      end if
      EnInSplit = EnFinSplit
    end do
    call cwtime(Ctime(2),Wtime(2))
    C_Dress = C_Dress+Ctime(2)-Ctime(1)
    W_Dress = W_Dress+Wtime(2)-Wtime(1)
    if (DBG) then
      !write(u6,*) 'CPU timing in UAA diag.: ',Ctime(2)-Ctime(1)
      write(u6,*) 'W. timing  in UAA diag: ',Wtime(2)-Wtime(1)
    end if

    if ((iterSplit == MxIterSplit) .and. (diffSplit > ThrSplit)) then
      if (ICIONLY == 0) then ! Hopefully the optimization will solve the convergence problem
        iErrSplit = 1
      else ! CIONLY case
        iErrSplit = 2
      end if
    end if
    if (iErrSplit /= 2) then
      if (DBG) then
        write(u6,*) 'After iteration :',iterSplit
        write(u6,*) 'In Root ',lrootSplit,' SplitCAS Energy :',EnInSplit
        write(u6,*) 'iDimBlockACNF from DiagOrd : ',iDimBlockACNF
        write(u6,*) 'iDimBlockA from DiagOrd : ',iDimBlockA
      end if
      !*****************************************************************
      ! coeffs c_m belonging to class (B)                              *
      !*****************************************************************
      ! TotSplitV will contain all the nConf coeff. for a single root.
      ! For a multi-root procedure the code should iterate over the roots.
      call mma_allocate(TotSplitV,nConf,label='totSplitVec')
      !TotSplitV(:) = Zero
      !call CmSplit(IPCSFtot,IPCNFtot,
      call cwtime(Ctime(1),Wtime(1))
      call get_Cm(IPCSFtot,IPCNFtot,nConf,NCNASM(STSYM),iDimBlockA,iDimBlockACNF,SplitV(:,lRootSplit),EnFinSplit,DTOC, &
                  DFTP,CONF,STSYM,HONE,ECORE,NAC,(NAEL+NBEL),NAEL,NBEL,TUVX,IPRINT,ExFac,IREOTS,FordSplit,TotSplitV)
      call cwtime(Ctime(2),Wtime(2))
      C_get_Cm = C_get_Cm+Ctime(2)-Ctime(1)
      W_get_Cm = W_get_Cm+Wtime(2)-Wtime(1)
      if (DBG) then
        write(u6,*) 'Get_Cm :'
        write(u6,*) 'CPU timing : ',Ctime(2)-Ctime(1)
        write(u6,*) 'W. timing  : ',Wtime(2)-Wtime(1)
      end if
      !*****************************************************************
      ! Normalization of the CI-Coefficients                           *
      !*****************************************************************
      SpliNor = ddot_(nConf,TotSplitV,1,TotSplitV,1)
      !write(u6,*) 'SpliNor',SpliNor
      TotSplitV(:) = TotSplitV(:)/sqrt(SpliNor)
      if (DBG) then
        call dVcPrt('Normalized...',' ',TotSplitV,nConf)
      end if
      !*****************************************************************
      ! SAVE CI VECTORS                                                *
      !*****************************************************************
      ! penso che in SplitCAS si debba usare sempre save_tmp_CI_vec.
      ! Ma nel dubbio lo copio anche con save_CI_vec:
      CIVEC(:) = Zero
      do j=1,nConf
        k = IPCSFtot(j)
        CIVEC(k) = TotSplitV(j)
      end do
      call Save_CI_vec(1,nConf,CIVEC,LuDavid)
      call Save_tmp_CI_vec(1,nConf,CIVEC,LuDavid)
      if (DBG) then
        write(String,'(A)') 'CI-diag in SplitCTL'
        call dVcPrt(String,' ',CIVEC,nConf)
      end if

      if (NAC == 0) then
        ENER(1,ITER) = EMY
      else
        !do jRoot=1,lRootSplit
        !  ENER(jRoot,ITER) = SplitE(jRoot)
        !  write(u6,*) 'ENER(jRoot,ITER)',ENER(jRoot,ITER)
        !end do
        ENER(lRootSplit,ITER) = SplitE(lRootSplit)
        if (DBG) then
          write(u6,*) 'lRootSplit :',lRootSplit
          write(u6,*) 'ITER :',ITER
          write(u6,*) 'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
        end if
      end if
    else
      write(u6,*) 'SplitCAS has to stop because it didn''t converge.'
    end if
    call xflush(u6)
    call mma_deallocate(TotSplitV)
    call mma_deallocate(Diag)
    call mma_deallocate(IPCSFtot)
    call mma_deallocate(IPCNFtot)
    !*******************************************************************
    ! CLEANUP AFTER SPLITCAS                                           *
    !*******************************************************************
    call mma_deallocate(AABlock)
    call mma_deallocate(DHAM)
    call mma_deallocate(SplitE)
    call mma_deallocate(SplitV)
    call mma_deallocate(HONE)
    call mma_deallocate(IPCNF)
    !call mma_deallocate(iSel)

  else
    !*******************************************************************
    ! Native Hamiltonian DIAGONALIZATION                               *
    !*******************************************************************
    ECORE = Zero
    MXSpli = iDimBlockA
    nAAblock = MXSpli*(MXSpli+1)/2
    call mma_allocate(HONE,NAC,NAC,label='HONE')
    call mma_allocate(iSel,MXSpli,label='iSel')
    call mma_allocate(IPCNF,NCNASM(STSYM),label='IPCNF')
    call mma_allocate(AABlock,nAAblock,label='AAblock')
    ! EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE
    call SQUARE(LW1,HONE,NAC,1,NAC)

    ! Calculate the AA Block of the Hamiltonian Matrix
    call mma_maxDBLE(MXXWS)
    call mma_allocate(Scr,MXXWS,label='EXHSCR')
    call PHPCSF(AABlock,iSel,IPCNF,MXSpli,DTOC,DFTP,CONF,STSYM,HONE,ECORE,NAC,Scr,NCNASM(STSYM),NAEL+NBEL,NAEL, &
                NBEL,iDimBlockA,iDimBlockACNF,CIVEC,TUVX,IPRINT,ExFac,IREOTS)
    call mma_deallocate(Scr)
    if (DBG) then
      call TRIPRT('AA block of the Hamiltonian Matrix',' ',AABlock,iDimBlockA)
    end if
    write(u6,*) '####################################################'
    write(u6,*) '# Dimension of AA block is equal to total nconf:   #'
    write(u6,*) '#            SplitCAS will not be used!            #'
    write(u6,*) '####################################################'
    iDimBlockTri = iDimBlockA*(iDimBlockA+1)/2
    call mma_allocate(DHAM,iDimBlockTri,label='DHAM')
    DHAM(:) = AABlock(:)

    call mma_allocate(SplitE,iDimBlockA,label='SplitE')
    call mma_allocate(SplitV,iDimBlockA,iDimBlockA,label='SplitV')
    !IPCSFtot(:) = iSel(:)
    call unitmat(SplitV,iDimBlockA)
    call NIdiag(DHAM,SplitV,iDimBlockA,iDimBlockA)
    call JACORD(DHAM,SplitV,iDimBlockA,iDimBlockA)
    do idx=1,iDimBlockA
      SplitE(idx) = DHAM(idx*(idx+1)/2)
    end do
    if (DBG) then
      call IVCPRT('Configurations included ',' ',iSel,iDimBlockA)
      call DVCPRT('Eigenval. of the explicit Hamiltonian',' ',SplitE,iDimBlockA)
      call RECPRT('Eigenvec. of the explicit Hamiltonian',' ',SplitV,iDimBlockA,iDimBlockA)
    end if
    !*******************************************************************
    ! SAVE CI VECTORS  after native Diagonalization                    *
    !*******************************************************************
    !write(u6,*) 'Root : ',lRootSplit
    !do i=1,lRootSplit
    CIVEC(:) = Zero
    do j=1,nConf
      k = iSel(j)
      CIVEC(k) = SplitV(j,lRootSplit)
    end do
    !call Save_tmp_CI_vec(i,lRootSplit,nConf,CIVEC,LuDavid)
    !call Save_tmp_CI_vec(1,lRootSplit,nConf,CIVEC,LuDavid)
    call Save_tmp_CI_vec(1,nConf,CIVEC,LuDavid)
    !if (IPRLEV == INSANE) then
    !  write(u6,'(A,I2)') 'Start vector of root',i
    !  write(u6,*) 'LuDavid',LuDavid
    !  write(String,'(A)') ' CI-coefficients in SplitCAS native'
    !  call dVcPrt(String,' ',CIVEC,nConf)
    !end if
    !end do
    if (NAC == 0) then
      ENER(1,ITER) = EMY
    else
      !do jRoot = 1,lRootSplit
      !  ENER(jRoot,ITER) = SplitE(jRoot)
      !  write(u6,*) 'ENER(jRoot,ITER)',ENER(jRoot,ITER)
      !end do
      ENER(lRootSplit,ITER) = SplitE(lRootSplit)
      if (DBG) then
        write(u6,*) 'ITER :',ITER
        write(u6,*) 'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
      end if
    end if
    !*******************************************************************
    ! CLEANUP AFTER native DIAGONALIZATION                             *
    !*******************************************************************
    call mma_deallocate(HONE)
    call mma_deallocate(IPCNF)
    call mma_deallocate(iSel)
    call mma_deallocate(AABlock)
    call mma_deallocate(DHAM)
    call mma_deallocate(SplitV)
    call mma_deallocate(SplitE)
  end if ! End of SplitCas/Native Hamiltonian diagonalization

else ! Do it IF there is CIRESTART
  iJOB = 0
  call f_Inquire('JOBOLD',Exists)
  if (Exists) iJOB = 1
  if (iJOB == 1) write(u6,*) ' Initial CI-vectors are read from JOBOLD'
  if (iJOB == 0) write(u6,*) ' Initial CI-vectors are read from JOBIPH'
  if (iJOB == 1) then
    if (JOBOLD <= 0) then
      JOBOLD = 20
      call DaName(JOBOLD,'JOBOLD')
    end if
  else
    JOBOLD = JOBIPH
  end if
  iDisk = 0
  call IDafile(JOBOLD,2,iToc,15,iDisk)
  iDisk = iToc(4)
  call mma_allocate(Tmp1,nConf,label='Scr1')
  call mma_allocate(Tmp2,nConf,label='Scr2')
  call mma_allocate(vkcnf,nactel,label='kcnf')
  !do i=1,lRootSplit
  call DDafile(JOBOLD,2,Tmp1,nConf,iDisk)
  call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp1,Tmp2,vkcnf)
  call Save_CI_vec(1,nConf,Tmp2,LuDavid)
  !write(u6,'(A,I2)') 'Start vector of root',i
  !if (DBG) then
  write(u6,*) 'LuDavid',LuDavid
  write(String,'(A)') '(CI coefficient in CIRST)'
  call dVcPrt(String,' ',Tmp2,nConf)
  !end if
  !end do
  call mma_deallocate(Tmp1)
  call mma_deallocate(Tmp2)
  call mma_deallocate(vkcnf)
  if (iJOB == 1) then
    if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
      call DaClos(JOBOLD)
      JOBOLD = -1
    else if (JOBOLD > 0) then
      JOBOLD = -1
    end if
  end if

end if ! End of do it IF there (is)/(isn't) CIRESTART

!-----------------------------------------------------------------------
!    CLEANUP AFTER SPLITCAS and DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------

! CIVEC: TEMPORARY CI VECTOR IN CSF BASIS

iDisk = IADR15(4)
!call Term_David(ICICH,ITERCI,lRootSplit,nConf,CIVEC,JOBIPH,LuDavid,iDisk)
!write(u6,*) 'ITERCI :',ITERCI
call Term_David(ICICH,ITERCI,1,nConf,CIVEC,JOBIPH,LuDavid,iDisk)
call mma_deallocate(CIVEC)

if (DBG) then
  call cwtime(CSplitTot(2),WSplitTot(2))
  write(u6,*) 'Total Time in SplitCAS:'
  write(u6,*) 'CPU timing : ',CSplitTot(2)-CSplitTot(1)
  write(u6,*) 'W. timing  : ',WSplitTot(2)-WSplitTot(1)
end if

return

end subroutine splitCTL
