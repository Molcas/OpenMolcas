************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      Subroutine PutRlx(D,DS,P,DAO,C)
      use spin_correlation, only: tRootGrad
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: CBLBM, ENER, ExFac, iBLBM, IPCMROOT, iPr,
     &                         iRLXRoot, iSymBB, ITER, lRoots, lSquare,
     &                         NACPAR, NACPR2, NewFock, nFint, nRoots,
     &                         NSXS, NTOT4, RlxGrd, iAdr15, ISTORP,
     &                         JBLBM
      use PrintLevel, only: DEBUG,USUAL
      use output_ras, only: LF,IPRLOC
      use general_data, only: NTOT1,NTOT2,NSYM,JOBIPH,NBAS
      use DWSol, only: DWSolv, DWSol_wgt, W_SOLV
      use rctfld_module, only: lRF

      Implicit None
      Real*8 D(*),DS(*),P(*),DAO(*),C(*)

      Character(LEN=8)  Fmt2
      Character(LEN=16), Parameter:: ROUTINE='PUTRLX  '
      Real*8 rdum(1)
      Integer i, iFinal, iPrLev, istmp, itmp, jDisk, jtmp, kDisk,
     &        left, NFSize, NZ
      Real*8 rTmp, wgt
      Real*8, External:: DNRM2_

      Real*8, Allocatable:: DA(:), DA_ave(:), DI(:), DSX(:), DS_ave(:),
     &                      DX(:), F(:), B(:), Q(:), FA(:), FI(:),
     &                      PUVX(:), TUVX(:), PX(:), WRK1(:), WRK2(:)
      Logical, external :: PCM_On

      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*
*     Read in corresponding density matrixes
*

      ! empty line before root resolved orbital gradient
      NZ=NTOT2
      kDisk = IADR15(3)  ! we need to save this address to reset later
      jDisk = IADR15(3)
      Call mma_allocate(DA,NZ,Label='DA')
      Call mma_allocate(DI,NZ,Label='DI')
      Call mma_allocate(DSX,NZ,Label='DSX')
      if (tRootGrad) then
        write(lf,*)
        do i = 1, LRoots
          call ddafile(JOBIPH, 2, D, NACPAR, jDisk)
          call ddafile(JOBIPH, 2, DS, NACPAR, jDisk)
          call ddafile(JOBIPH, 2, P, NACPR2, jDisk)
          call ddafile(JOBIPH, 0, rdum, NACPR2, jDisk)
*
*         Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
*
          Call Get_D1I_RASSCF(C,DI)
*
          Call mma_allocate(DX,NACPAR,Label='DX')
*
          CALL DCOPY_(NACPAR,DS,1,DX,1)
          Call DBLOCK(DX)
          CALL Get_D1A_RASSCF(C,DX,DSX)
*
          CALL DCOPY_(NACPAR,D,1,DX,1)
          Call DBLOCK(DX)
          CALL Get_D1A_RASSCF(C,DX,DA)
*
*         Construct the Fock matrix used for the connection term.
*
          NFSIZE=MAX(NACPAR,NTOT4)
          Call mma_allocate(F,NFSIZE,Label='F')
          Call mma_allocate(B,NZ,Label='B')
          Call mma_allocate(Q,NZ,Label='Q')
          Call mma_allocate(FA,NZ,Label='FA')
          Call mma_allocate(FI,NZ,Label='FI')
          Call mma_allocate(PUVX,NFINT,Label='PUVX')
          Call mma_allocate(tuvx,nacpr2,Label='tuvx')
          IPR=0
          IF(IPRLEV.EQ.3) IPR=1
          IF(IPRLEV.EQ.4) IPR=5
          IF(IPRLEV.EQ.5) IPR=10
          PUVX(:)=0.0D0
          Call TraCtl2(C,PUVX,TUVX,DI,FI,DA,FA,ipr,lsquare,ExFac)
          CALL SGFCIN(C,F,FI,DI,DA,DSX)
          call dcopy_(ntot4,[0.0d0],0,F,1)
          call dcopy_(ntot4,[0.0d0],0,B,1)
*
*         Prevent FMAT from changing Active fock matrix
*
          iTmp=newfock
          newFock=-99999
*
          Call Fmat(C,PUVX,DX,DA,FI,FA)
          newFock=itmp
          IFINAL=1
          rTmp=   CBLBM
          itmp= IBLBM
          jtmp= JBLBM
          istmp= ISYMBB

          IF(ISTORP(NSYM+1).GT.0) THEN
              CALL mma_allocate(PX,ISTORP(NSYM+1),Label='PX')
              CALL PMAT_RASSCF(P,PX)
          Else
              CALL mma_allocate(PX,1,Label='PX')
          END IF
*
          Call FOCK(F,B,FI,FA,DX,PX,Q,PUVX,IFINAL,C)
          CBLBM=rtmp
          IBLBM=itmp
          JBLBM=jtmp
          ISYMBB=istmp
          CALL mma_deallocate(PX)

          write(6,'(6x,a,i3,5x,f12.10)')
     &        "Norm of electronic gradient for root ", i,
     &        DNRM2_(NSXS,B,1)
*
          Call mma_deallocate(DX)
          Call mma_deallocate(tuvx)
          Call mma_deallocate(puvx)
          Call mma_deallocate(FI)
          Call mma_deallocate(FA)
          Call mma_deallocate(Q)
          Call mma_deallocate(B)
          Call mma_deallocate(F)
        end do
      end if

      Do i=1,iRlxRoot-1
        Call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
      End Do
      Call DDaFile(JOBIPH,2,D,NACPAR,kDisk)
      Call DDaFile(JOBIPH,2,DS,NACPAR,kDisk)
      Call DDaFile(JOBIPH,2,P,NACPR2,kDisk)
      Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
*
*     Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
      Call Get_D1I_RASSCF(C,DI)
*
      Call mma_allocate(DX,NACPAR,Label='DX')
*
      CALL DCOPY_(NACPAR,DS,1,DX,1)
      Call DBLOCK(DX)
      CALL Get_D1A_RASSCF(C,DX,DSX)
*
      CALL DCOPY_(NACPAR,D,1,DX,1)
      Call DBLOCK(DX)
      CALL Get_D1A_RASSCF(C,DX,DA)
*
*     Construct the Fock matrix used for the connection term.
*
      NFSIZE=MAX(NACPAR,NTOT4)
      Call mma_allocate(F,NFSIZE,Label='F')
      Call mma_allocate(B,NZ,Label='B')
      Call mma_allocate(Q,NZ,Label='Q')
      Call mma_allocate(FA,NZ,Label='FA')
      Call mma_allocate(FI,NZ,Label='FI')
      Call mma_allocate(PUVX,NFINT,Label='PUVX')
      Call mma_allocate(tuvx,nacpr2,Label='tuvx')
      IPR=0
      IF(IPRLEV.EQ.3) IPR=1
      IF(IPRLEV.EQ.4) IPR=5
      IF(IPRLEV.EQ.5) IPR=10
      PUVX(:)=0.0D0
      Call TraCtl2(C,PUVX,TUVX,DI,FI,DA,FA,ipr,lsquare,ExFac)
*
      ! DA constructed above is the density for the RLXROOT state.
      ! We should reconstruct the AO density matrix to prevent the
      ! mismatch of states for reaction field and geometry optimization
      ! It seems that IPCMROOT may not be defined sometimes
      if (lRF .and. PCM_On() .and. (IPCMROOT /= iRLXROOT .or.
     *    IPCMROOT <= 0 .or. DWSolv%DWZeta /= 0.0d+00)) then
        !! Polarize PCM etc with weighted density
        Call mma_allocate(DA_ave,MAX(NACPAR,NZ),Label='DA_ave')
        Call mma_allocate(DS_ave,MAX(NACPAR,NZ),Label='DS_ave')
        Call mma_allocate(WRK1,MAX(NACPAR,NZ),Label='WRK1')
        Call mma_allocate(WRK2,MAX(NACPAR,NZ),Label='WRK2')
*
        DA_ave(:) = 0.0d+00
        DS_ave(:) = 0.0d+00
*
        call DWSol_wgt(2,ENER(:,ITER))
        kDisk = IADR15(3)
        Do i=1,nRoots
          wgt = W_SOLV(i)
          Call DDaFile(JOBIPH,2,WRK1,NACPAR,kDisk)
          Call DDaFile(JOBIPH,2,WRK2,NACPAR,kDisk)
          Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
          Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
          if (wgt < 1.0d-09) cycle
          call daxpy_(NACPAR,wgt,WRK1,1,DA_ave,1)
          call daxpy_(NACPAR,wgt,WRK2,1,DS_ave,1)
        End Do
*
*       Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
        CALL DCOPY_(NACPAR,DS_ave,1,DX,1)
        call dcopy_(NZ,DSX,1,DS_ave,1)
        Call DBLOCK(DX)
        CALL Get_D1A_RASSCF(C,DX,DSX)
*
        CALL DCOPY_(NACPAR,DA_ave,1,DX,1)
        call dcopy_(NZ,DA,1,DA_ave,1)
        Call DBLOCK(DX)
        CALL Get_D1A_RASSCF(C,DX,DA)
*
        Call mma_deallocate(WRK1)
        Call mma_deallocate(WRK2)
        IF(IPRLEV.GE.USUAL .AND. DWSolv%DWZeta > 0.0d+00) THEN
          left=6
          Write(LF,*)
          Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
          Write(LF,Fmt2//'A)')
     &      'Dynamically weighted solvation has been employed'
          Write(LF,Fmt2//'A,(T45,10F6.3))')
     &      'Final weights for the reaction field:',
     &                                (W_SOLV(i),i=1,nRoots)
        END IF
      end if
*
      CALL SGFCIN(C,F,FI,DI,DA,DSX)
      call dcopy_(ntot4,[0.0d0],0,F,1)
      call dcopy_(ntot4,[0.0d0],0,B,1)
*
*     Prevent FMAT from changing Active fock matrix
*
      iTmp=newfock
      newFock=-99999
*
      if (lRF .and. PCM_On() .and. (IPCMROOT /= iRLXROOT .or.
     *    IPCMROOT <= 0 .or. DWSolv%DWZeta /= 0.0d+00)) then
        ! The rest of the RASSCF program uses state-specific density
        ! (note that, it is iRlxRoot!), so restore the one constructed
        ! above, before TraCtl2
        Call dcopy_(NZ,DA_ave,1,DA,1)
        Call dcopy_(NZ,DS_ave,1,DSX,1)
        Call mma_deallocate(DA_ave)
        Call mma_deallocate(DS_ave)
*
        Call mma_allocate(WRK1,nTot1,Label='WRK1')
        Call mma_allocate(WRK2,nTot1,Label='WRK2')
        Call Fold(nSym,nBas,DI,WRK1)
        Call Fold(nSym,nBas,DA,WRK2)
        Call Daxpy_(nTot1,1.0d+00,WRK1,1,WRK2,1)
        Call Put_dArray('D1ao',WRK2,nTot1)
        Call Fold(nSym,nBas,DSX,WRK1)
        Call Put_dArray('D1sao',WRK1,nTot1)
        Call mma_deallocate(WRK1)
        Call mma_deallocate(WRK2)
      end if
*
      Call Fmat(C,PUVX,DX,DA,FI,FA)
      newFock=itmp
      IFINAL=1
      rTmp=   CBLBM
      itmp= IBLBM
      jtmp= JBLBM
      istmp= ISYMBB

      IF(ISTORP(NSYM+1).GT.0) THEN
          CALL mma_allocate(PX,ISTORP(NSYM+1),Label='PX')
          CALL PMAT_RASSCF(P,PX)
      Else
          CALL mma_allocate(PX,1,Label='PX')
      END IF
*
      Call FOCK(F,B,FI,FA,DX,PX,Q,PUVX,IFINAL,C)
      CBLBM=rtmp
      IBLBM=itmp
      JBLBM=jtmp
      ISYMBB=istmp
      CALL mma_deallocate(PX)

      RlxGrd=DNRM2_(NSXS,B,1)
*
      Call mma_deallocate(DX)
      Call mma_deallocate(tuvx)
      Call mma_deallocate(puvx)
      Call mma_deallocate(FI)
      Call mma_deallocate(FA)
      Call mma_deallocate(Q)
      Call mma_deallocate(B)
      Call mma_deallocate(F)
*
* Add up one electron densities
*
      call daxpy_(nZ,1.0d0,DA,1,DI,1)
      Call Fold(nSym,nBas,DI,DAO)

      Call mma_deallocate(DSX)
      Call mma_deallocate(DA)
      Call mma_deallocate(DI)

      end Subroutine PutRlx
