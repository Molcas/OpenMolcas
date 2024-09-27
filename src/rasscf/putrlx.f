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
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "output_ras.fh"
      Character(LEN=16), Parameter:: ROUTINE='PUTRLX  '
      Real*8 D(*),DS(*),P(*),DAO(*),C(*)
      Real*8 rdum(1)

      Real*8, Allocatable:: DA(:), DI(:), DSX(:), DX(:), F(:), B(:),
     &                      Q(:), FA(:), FI(:), PUVX(:), TUVX(:), PX(:)

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
      CALL SGFCIN(C,F,FI,DI,DA,DSX)
      call dcopy_(ntot4,[0.0d0],0,F,1)
      call dcopy_(ntot4,[0.0d0],0,B,1)
*
*     Prevent FMAT from changing Active fock matrix
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
