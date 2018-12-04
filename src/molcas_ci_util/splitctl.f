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
* Copyright (C) 2009, Giovanni Li Manni                                *
*               2009, Francesco Aquilante                              *
************************************************************************

      Subroutine splitCTL(LW1, TUVX, IFINAL,iErrSplit)
************************************************************************
*     CI Hamiltonian Matrix elements reader                            *
*     calling arguments:                                               *
*     LW1     : Memory pointer to active Fock matrix                   *
*               array of real*8                                        *
*     TUVX    : array of real*8                                        *
*               two-electron integrals (tu!vx)                         *
*     IFINAL  : integer                                                *
*               termination flag                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by:                                                      *
*     G. Li Manni (GLMJ) and F. Aquilante                              *
*     University of Geneva, Switzerland, 2009                          *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)
      Character*80 String
      Logical DBG, Exist
      Dimension LW1(*), TUVX(*)

#include "rasdim.fh"
#include "rasscf.fh"
#include "splitcas.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter(ROUTINE='splitCTL')
#include "csfbas.fh"
#include "spinfo.fh"
#include "strnum.fh"
#include "timers.fh"

      Call qEnter('SplitCTL')

      IPRLEV=IPRLOC(3)

      DBG=.false.
      DBG=IPRLEV.GE.DEBUG

      iErrSplit = 0
* -------------------------------------------------------------------- C
* -- INITIALIZE THE DAVIDSON DIAGONALIZATION
* -------------------------------------------------------------------- C
*
*1) In posizione 1 nella seguente routine ho messo 'lRootSplit' invece
*   di 'lRoots', in quanto 'lRoots' viene specificato nell'input del
*   calcolo con la keyword CIROOT (secondo numero intero scritto nella
*   linea sotto CIROOT) ma sottomettendo un calcolo SplitCAS non andiamo
*   ad usare tale keyword (almeno non per il momento!).
*
*2) In posizione 4 nella seguente routine ho messo 'nConf' invece di 'nSel'
*   rispetto alla routine originale in 'davctl'.
*
* Qui nConf rappresenta il numero totale di CSFs.

      call cwtime(CSplitTot1,WSplitTot1)
*     Call Ini_David(lRootSplit,nConf,nDet,nconf,nAc,LuDavid)
      Call Ini_David(1,nConf,nDet,nconf,nAc,LuDavid)
* -------------------------------------------------------------------- C
* --  COMPUTE THE DIAGONAL ELEMENTS OF THE COMPLETE HAMILTONIAN
* -------------------------------------------------------------------- C
*           LW4: TEMPORARY CI VECTOR IN CSF BASIS
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)

        goto 29555
        If(IFINAL.eq.2) then
*       ^to avoid last diagonalization
          Call dCopy_(nConf,0.0d0,0,Work(LW4),1)
*          Call Load_tmp_CI_vec(1,1,nConf,Work(LW4),LuDavid)
           Call Load_CI_vec(1,nConf,Work(LW4),LuDavid)
*          Call dDaFile(JOBIPH,2,Work(LW4),nConf,LuDavid)
          if (DBG) then
            write(6,*) 'LuDavid', LuDavid
            Write (String,'(A)') 'Final=2 : CI-coeff in SplitCAS'
            Call dVcPrt(String,' ',Work(LW4),nConf)
          end if

          If ( NAC.eq.0 ) then
            ENER(1,ITER)=EMY
          Else
*              ENER(lRootSplit,ITER) = ENER(lRootSplit,ITER-1)
              if (DBG) then
                write(6,*)'lRootSplit :',lRootSplit
                write(6,*)'ITER :',ITER
                write(6,*)'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
              end if
          End if
          CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)
          Return
        End if
29555 continue

      IF(NAC. GT. 0)
     &  CALL CIDIA_CI_UTIL(NAC,NCONF,LSYM,WORK(LW4),LW1,TUVX,
     &                          LUDAVID)
************************************************************************
* iCaseSplit = 1  : there is NOT CI-RESTART.
* iCaseSplit = 2  : there is CIRESTART. The code will read the CI
*                   coefficients from the JOBOLD file.
************************************************************************
      iCaseSplit = 1
      if (ICIRST.ne.0) iCaseSplit = 2

      IF (iCaseSplit.eq.1) THEN
*     ^ There is NO CIRST
      if (iDimBlockA.ne.nconf) then
*     ^ AA Block is smaller than the full Hamiltonian matrix, then
*       SplitCAS will be performed.

        CALL GETMEM('HONE','ALLO','REAL',LOCONE,NAC**2)
* EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE
        CALL TRIEXP(LW1,Work(LOCONE),NAC)
        CALL GETMEM('IREOTS','ALLO','INTEGER',IREOTS,NAC)
        CALL GET_IREOTS(IWORK(IREOTS),NAC)
C        CALL GETMEM('IPCNF','ALLO','INTE',LG1,NCNASM(LSYM))

        call cwtime(C_ABlockDim_sel1,W_ABlockDim_sel1)
        ECORE=0.0D0
        if ( NumSplit ) then
************************************************************************
* In such case the splitting is done according to a NUMERICAL          *
* selection.                                                           *
* In ipCSFSplit the following quantities are calculated:               *
*    1) The index array for the CSFs and CNFs                          *
*    2) The final values of iDimBlockA and iDimBlockACNF               *
*    3) The diagonal array of the Hamiltonian matrix                   *
************************************************************************
        CALL GETMEM('DiagCSF','ALLO','REAL',ipDiag,nConf)
        Call GetMem('IPCSFtot','Allo','Integer',ipCSFtot,nconf)
        Call GetMem('IPCNFtot','Allo','Integer',ipCNFtot,NCNASM(LSYM))
        CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
        CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)
        MXSpli = iDimBlockA
*        nAAblock = MXSpli*(MXSpli+1)/2
        CALL ipCSFSplit(Work(ipDiag),iWork(ipCSFtot),iWork(ipCNFtot),
     &                  nConf,MXSpli,
     &                  Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),
     &                  LSYM,Work(LOCONE),ECORE,NAC,
     &                  Work(LW2),NCNASM(LSYM),(NAEL+NBEL),NAEL,NBEL,
     &                  Work(LW4),TUVX,
     &                  IPRINT,ExFac,IWORK(IREOTS))
*        CALL DVCPRT('Diagonal elements of Hamilt. matrix in CSF',' ',
*     &              Work(ipDiag),nConf)
        CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
*        CALL GETMEM('DiagCSF','FREE','REAL',ipDiag,nConf)
*        nAAblock = iDimBlockA*(iDimBlockA+1)/2
        end if

        if ( EnerSplit .or. PerSplit ) then
************************************************************************
*    COMPUTE the index array for THE DIAGONAL ELEMENTS OF THE          *
*    COMPLETE HAMILTONIAN in energetic order.                          *
************************************************************************
        CALL GETMEM('DiagCSF','ALLO','REAL',ipDiag,nConf)
        CALL GETMEM('DiagCNF','ALLO','REAL',ipDiagCNF,NCNASM(LSYM))
        CALL GetMem('ipCSFtot','Allo','Integer',ipCSFtot,nConf)
        CALL GETMEM('ipCNFtot','ALLO','Integer',ipCNFtot,NCNASM(LSYM))
        CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
        CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)
* 'GapSpli' comes from the input in eV
* 'condition' goes to DiagOrd in Hartree if EnerSplit
* 'condition' goes to DiagOrd as a percentage if PerSplit
        if (EnerSplit ) condition = gapSpli/27.2114
        if ( PerSplit ) condition = percSpli
        call DiagOrd(Work(ipDiag),Work(ipDiagCNF),
     &             iWork(ipCSFtot),iWork(ipCNFtot),nConf,condition,ITER,
     &              Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),
     &              LSYM,Work(LOCONE),ECORE,NAC,
     &              Work(LW2),NCNASM(LSYM),(NAEL+NBEL),NAEL,NBEL,
     &              TUVX,IPRINT,ExFac,IWORK(IREOTS))
        if (DBG) then
          CALL DVCPRT('Diagonal elements of Hamilt. matrix in CSF',' ',
     &                  Work(ipDiag),nConf)
          CALL DVCPRT('Diagonal elements of Hamilt. matrix in CNF',' ',
     &                  Work(ipDiagCNF),NCNASM(LSYM))
          CALL IVCPRT('Index Array in CSF',' ',
     &                  iWork(ipCSFtot),nConf)
          CALL IVCPRT('Index Array in CNF',' ',
     &                  iWork(ipCNFtot),NCNASM(LSYM))
          write(6,*) 'iDimBlockACNF from DiagOrd : ', iDimBlockACNF
          write(6,*) 'iDimBlockA from DiagOrd : ', iDimBlockA
        end if
        call xflush(6)
        CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
        CALL GETMEM('DiagCNF','free','REAL',ipDiagCNF,NCNASM(LSYM))
*        CALL GETMEM('DiagCSF','free','REAL',ipDiag,nConf)
*        CALL GETMEM('indOrdCNF','free','Integer',ipCNFtot,NCNASM(LSYM))
        end if

        if (DBG) then
         call cwtime(C_ABlockDim_sel2,W_ABlockDim_sel2)
         write(6,*) 'Time needed to select CSFs in A-Block:'
         write(6,*) 'CPU timing : ', C_ABlockDim_sel2 - C_ABlockDim_sel1
         write(6,*) 'W. timing  : ', W_ABlockDim_sel2 - W_ABlockDim_sel1
        end if
************************************************************************
*  Let's start with SPLITCAS code                                      *
************************************************************************
*        call cwtime(C_iterative1,W_iterative1)
        call chksplit()
************************************************************************
      if (DBG) then
        write(6,*) 'iDimBlockA after selection : ', iDimBlockA
        write(6,*) 'iDimBlockACNF after selection : ', iDimBlockACNF
      end if
* calculate the dressed HAMILTONIAN matrix (STEP 2)
        iDimBlockTri= iDimBlockA*(iDimBlockA+1)/2
*        iBCSF = nConf-iDimBlockA
*        iBCNF = NCNASM(LSYM)-NPCNF
*        iBMAX = MAX(iBCSF,iDimBlockA)
        CALL GETMEM('AAblock','ALLO','REAL',ipAABlock,iDimBlockTri)
        CALL GETMEM('DHAM','ALLO','REAL',ipDHAM,iDimBlockTri)
*        CALL GETMEM('BVEC','ALLO','REAL',ipBVEC,iBMAX)
        call Fzero(Work(ipAABlock),  iDimBlockTri)
        if(iter.eq.1) then
          EnInSplit=Work(ipDiag+lrootSplit-1)
        end if
        if (DBG) then
          write(6,*) 'Initial Energy in SplitCAS : ', EnInSplit
        end if
        diffSplit = ThrSplit + 1.0d0

        iterSplit = 0
        call cwtime(C_Dress_1,W_Dress_1)
        call getmem('SplitE','Allo','Real',ipSplitE,iDimBlockA)
        call getmem('SplitV','Allo','Real',ipSplitV,
     &              iDimBlockA*iDimBlockA)
        DO WHILE (diffSplit.gt.ThrSplit.and.iterSplit.lt.MxIterSplit)
          iterSplit = iterSplit + 1
          if (DBG) then
            write(6,*) '*************** Iteration SplitCAS =', iterSplit
          end if
C          call Compute_Umn(Work(ipBVEC),NPCNF,NCNASM(LSYM),
C     &                     EnInSplit,NPCNF+1,1,Work(ipDHAM))
C         CALL SPLITCSF(Work(ipAABlock),EnInSplit,Work(ipDHAM),
          CALL get_Umn(Work(ipAABlock),EnInSplit,Work(ipDHAM),
     &                  iWork(ipCSFtot),iWork(ipCNFtot),nconf,
     &                  Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),
     &                  LSYM,Work(LOCONE),ECORE,NAC,
     &                  NCNASM(LSYM),(NAEL+NBEL),NAEL,NBEL,
     &                  iDimBlockA,iDimBlockACNF,Work(LW4),TUVX,
     &                  iterSplit,ITER,
     &                  IPRINT,ExFac,IWORK(IREOTS))
          call xflush(6)
          if (DBG) then
            call TRIPRT('AA block of the Hamiltonian Matrix',' ',
     &                 Work(ipAABLOCK),iDimBlockA)
            call TRIPRT('Dressed AA block Hamiltonian Matrix',' ',
     &                  Work(ipDHAM),iDimBlockA)
            call xflush(6)
          end if
*         CALL GETMEM('BVEC','FREE','REAL',ipBVEC,iBMAX)
************************************************************************
*      Dressed Hamiltonian DIAGONALIZATION                                *
************************************************************************
          call dcopy_(iDimBlockA*iDimBlockA,0.0d0,0,Work(ipSplitV),1)
          do i=1,iDimBlockA
            ii=i+iDimBlockA*(i-1)
            Work(ipSplitV+ii-1)=1.0D00
          end do
          Call NIdiag(Work(ipDHAM),Work(ipSplitV),
     &                iDimBlockA,iDimBlockA,0)
          Call JACORD(Work(ipDHAM),Work(ipSplitV),
     &                iDimBlockA,iDimBlockA)
          do idx = 1, iDimBlockA
            Work(ipSplitE+idx-1) = Work(ipDHAM-1+idx*(idx+1)/2)
          end do
          EnFinSplit= Work(ipSplitE+lrootSplit-1)
          if (DBG) then
            CALL IVCPRT('CSFs included ',' ',iWork(ipCSFtot),iDimBlockA)
            CALL RECPRT('Eigenvec. of dressed Hamiltonian',' ',
     &                  Work(ipSplitV),iDimBlockA,iDimBlockA)
            CALL DVCPRT('Eigenval. of dressed Hamiltonian',' ',
     &                  Work(ipSplitE),iDimBlockA)
          end if
          diffSplit = abs(EnFinSplit-EnInSplit)
          if (DBG) write(6,*)'Energy diff in SplitCAS :', diffSplit
          if (DBG) then
            if (iterSplit.eq.1)
     &          write(6,*) 'IterSplit   Final Energy in SplitCAS '
            write(6,*) IterSplit, EnFinSplit
          end if
          EnInSplit = EnFinSplit
        END DO
          call cwtime(C_Dress_2,W_Dress_2)
          C_Dress_3 = C_Dress_3 + C_Dress_2 - C_Dress_1
          W_Dress_3 = W_Dress_3 + W_Dress_2 - W_Dress_1
*          write(6,*) 'CPU timing in UAA diag.: ', C_Dress_2 - C_Dress_1
*        if (DBG)
*     &    write(6,*) 'CPU timing in UAA diag.: ', C_Dress_2 - C_Dress_1
        if (DBG)
     &     write(6,*) 'W. timing  in UAA diag: ', W_Dress_2 - W_Dress_1

        if(iterSplit.eq.MxIterSplit.and.diffSplit.gt.ThrSplit) then
          if(ICIONLY.eq.0) then
*         ^Hopefully the optimization will solve the convergence problem
            iErrSplit = 1
          else
*         ^ CIONLY case
            iErrSplit = 2
          end if
        end if
        if(iErrSplit.ne.2) then
          if (DBG) then
            write(6,*)'After iteration :', iterSplit
            write(6,*)'In Root ',lrootSplit,' SplitCAS Energy :',
     &              EnInSplit
            write(6,*) 'iDimBlockACNF from DiagOrd : ', iDimBlockACNF
            write(6,*) 'iDimBlockA from DiagOrd : ', iDimBlockA
          end if
************************************************************************
* coeffs c_m belonging to class (B)                                    *
************************************************************************
*       iBBlockDim = nconf - iDimBlockA
*      call getmem('HmmDen','allo','real',ipHmmD,iBBlockDim)
*      call getmem('UmnCnold','allo','real',ipUCOld,iBBlockDim)
*      call getmem('UmnCnnew','allo','real',ipUCnew,iBBlockDim)
* Work(ipTotSplitV) will contain all the nConf coeff. for a single root.
* For a multi-root procedure the code should iterate over the roots.
          call getmem('totSplitVec','allo','real',ipTotSplitV,nConf)
*         call Fzero(Work(ipTotSplitV),nConf)
*         call CmSplit(iWork(ipCSFtot),iWork(ipCNFtot),
          call cwtime(C_get_Cm1,W_get_Cm1)
          call get_Cm (iWork(ipCSFtot),iWork(ipCNFtot),
     &                 nConf,NCNASM(LSYM),
     &                 iDimBlockA,iDimBlockACNF,
     &                 Work(ipSplitV+(lRootSplit-1)*iDimBlockA),
     &                 lrootSplit,EnFinSplit,
     &                 Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),
     &                 LSYM,Work(LOCONE),ECORE,
     &                 NAC,(NAEL+NBEL),NAEL,NBEL,
     &                 Work(LW4),TUVX,
     &                 IPRINT,ExFac,IWORK(IREOTS),
     &                 FordSplit,
     &                 Work(ipTotSplitV))
          call cwtime(C_get_Cm2,W_get_Cm2)
          C_get_Cm3 = C_get_Cm3 + C_get_Cm2 - C_get_Cm1
          W_get_Cm3 = W_get_Cm3 + W_get_Cm2 - W_get_Cm1
          if (DBG) then
            write(6,*) 'Get_Cm :'
            write(6,*) 'CPU timing : ', C_get_Cm2 - C_get_Cm1
            write(6,*) 'W. timing  : ', W_get_Cm2 - W_get_Cm1
          end if
*      call getmem('UmnCnnew','Free','real',ipUCnew,iBBlockDim)
*      call getmem('UmnCnold','Free','real',ipUCOld,iBBlockDim)
*      call getmem('HmmDen','Free','real',ipHmmD,iBBlockDim)
************************************************************************
*       Normalization of the CI-Coefficients                           *
************************************************************************
         SpliNor  = ddot_(nConf,Work(ipTotSplitV),1,Work(ipTotSplitV),1)
*       write(6,*)'SpliNor', SpliNor
         SqSpliNor= 1.0d0/sqrt(SpliNor)
         call dscal_(nConf,SqSpliNor,Work(ipTotSplitV),1)
         if (DBG) then
           Call dVcPrt('Normalized...',' ',Work(ipTotSplitV),nConf)
         end if
************************************************************************
************************************************************************
*      SAVE CI VECTORS                                                 *
************************************************************************
* penso che in SplitCAS si debba usare sempre save_tmp_CI_vec.
* Ma nel dubbio lo copio anche con save_CI_vec:
            Call dCopy_(nConf,0.0d0,0,Work(LW4),1)
            do j =1,nConf
              k= iWork(ipCSFtot+j-1)
              Work(LW4+k-1) = Work(ipTotSplitV+j-1)
            end do
            Call   Save_CI_vec(1,nConf,Work(LW4),LuDavid)
            Call Save_tmp_CI_vec(1,nConf,Work(LW4),LuDavid)
            if (DBG) then
              Write (String,'(A)') 'CI-diag in SplitCTL'
              Call dVcPrt(String,' ',Work(LW4),nConf)
            end if

          If ( NAC.eq.0 ) then
            ENER(1,ITER)=EMY
          Else
*           Do jRoot = 1,lRootSplit
*             ENER(jRoot,ITER) = Work(ipSplitE+jRoot-1)
              ENER(lRootSplit,ITER) = Work(ipSplitE+lRootSplit-1)
*             write(6,*)'ENER(jRoot,ITER)',ENER(jRoot,ITER)
              if (DBG) then
                write(6,*)'lRootSplit :',lRootSplit
                write(6,*)'ITER :',ITER
                write(6,*)'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
              end if
*           End Do
          End if
        else
         write(6,*)'SplitCAS has to stop because it didn''t converge.'
        end if
        call xflush(6)
        call getmem('totSplitVec','free','real',ipTotSplitV,nConf)
*        CALL GetMem('indOrdCSF','free','Integer',ipCSFtot,nConf)
        Call GetMem('IPCSFtot','FREE','Integer',ipCSFtot,nconf)
        Call GetMem('IPCNFtot','FREE','Integer',ipCNFtot,NCNASM(LSYM))
        CALL GetMem('DiagCSF','FREE','REAL',ipDiag,nConf)
************************************************************************
*            CLEANUP AFTER SPLITCAS                                    *
************************************************************************
        call getmem('SplitE','FREE','Real',ipSplitE,iDimBlockA)
        call getmem('SplitV','FREE','Real',ipSplitV,
     &                                     iDimBlockA*iDimBlockA)
        CALL GETMEM('DHAM','FREE','REAL',ipDHAM,iDimBlockTri)
        CALL GETMEM('AAblock','FREE','REAL',ipAABlock,iDimBlockTri)
        CALL GETMEM('IREOTS','FREE','INTEGER',IREOTS,NAC)
        CALL GETMEM('HONE','FREE','REAL',LOCONE,NAC**2)
C        CALL GETMEM('IPCNF','FREE','INTE',LG1,NCNASM(LSYM))
*        Call GetMem('iSel','Free','Integer',lSel,MxSpli)

      else
************************************************************************
* Native Hamiltonian DIAGONALIZATION                                   *
************************************************************************
      ECORE=0.0D0
      MXSpli=iDimBlockA
      nAAblock = MXSpli*(MXSpli+1)/2
      CALL GETMEM('HONE','ALLO','REAL',LOCONE,NAC**2)
      Call GetMem('iSel','Allo','Integer',lSel,MxSpli)
      CALL GETMEM('IPCNF','ALLO','INTE',LG1,NCNASM(LSYM))
      CALL GETMEM('AAblock','ALLO','REAL',ipAABlock,nAAblock)
* EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE
      CALL TRIEXP(LW1,Work(LOCONE),NAC)

* Calculate the AA Block of the Hamiltonian Matrix
      CALL GETMEM('IREOTS','ALLO','INTEGER',IREOTS,NAC)
      CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
      CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)
      CALL GET_IREOTS(IWORK(IREOTS),NAC)
      CALL PHPCSF(Work(ipAABlock),iWork(lSel),iWork(LG1),MXSpli,
     &              Work(KDTOC),iWork(KDFTP),iWork(KICONF(1)),
     &              LSYM,Work(LOCONE),ECORE,NAC,
     &              Work(LW2),NCNASM(LSYM),(NAEL+NBEL),NAEL,NBEL,
     &              iDimBlockA,iDimBlockACNF,Work(LW4),TUVX,
     &              IPRINT,ExFac,IWORK(IREOTS))
      CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
      if (DBG) then
        call TRIPRT('AA block of the Hamiltonian Matrix',' ',
     &              Work(ipAABLOCK),iDimBlockA)
      End if
        write(6,*)'####################################################'
        write(6,*)'# Dimension of AA block is equal to total nconf:   #'
        write(6,*)'#            SplitCAS will not be used!            #'
        write(6,*)'####################################################'
        iDimBlockTri= iDimBlockA*(iDimBlockA+1)/2
        CALL GETMEM('DHAM','ALLO','REAL',ipDHAM,iDimBlockTri)
        call dcopy_(iDimBlockTri,Work(ipAABlock),1,Work(ipDHAM),1)

        call getmem('SplitE','Allo','Real',ipSplitE,iDimBlockA)
        call getmem('SplitV','Allo','Real',ipSplitV,
     &              iDimBlockA*iDimBlockA)
        call dcopy_(iDimBlockA*iDimBlockA,0.0d0,0,Work(ipSplitV),1)
*        call icopy(iDimBlockA,iWork(lSel),1,iWork(ipCSFtot),1)
        do i=1,iDimBlockA
          ii=i+iDimBlockA*(i-1)
          Work(ipSplitV+ii-1)=1.0D00
        end do
        Call NIdiag(Work(ipDHAM),Work(ipSplitV),iDimBlockA,iDimBlockA,0)
        Call JACORD(Work(ipDHAM),Work(ipSplitV),iDimBlockA,iDimBlockA)
        do idx = 1, iDimBlockA
          Work(ipSplitE+idx-1) = Work(ipDHAM-1+idx*(idx+1)/2)
        end do
        if (DBG) then
          CALL IVCPRT('Configurations included ',' ',
     &               iWork(lSel),iDimBlockA)
          CALL DVCPRT('Eigenval. of the explicit Hamiltonian',' ',
     &               Work(ipSplitE),iDimBlockA)
          CALL RECPRT('Eigenvec. of the explicit Hamiltonian',' ',
     &               Work(ipSplitV),iDimBlockA,iDimBlockA)
        End if
************************************************************************
*      SAVE CI VECTORS  after native Diagonalization                   *
************************************************************************
*        Write (6,*) 'Root : ',lRootSplit
*       Do i = 1,lRootSplit
          Call dCopy_(nConf,0.0d0,0,Work(LW4),1)
          do j =1,nConf
            k= iWork(lSel+j-1)
            Work(LW4+k-1) = Work(ipSplitV+(lRootSplit-1)*nconf+j-1)
          end do
*         Call Save_tmp_CI_vec(i,lRootSplit,nConf,Work(LW4),LuDavid)
*         Call Save_tmp_CI_vec(1,lRootSplit,nConf,Work(LW4),LuDavid)
          Call Save_tmp_CI_vec(1,nConf,Work(LW4),LuDavid)
*         If ( IPRLEV.eq. INSANE ) then
*           Write (6,'(A,I2)') 'Start vector of root',i
*            write(6,*)'LuDavid',LuDavid
*            Write (String,'(A)') ' CI-coefficients in SplitCAS native'
*            Call dVcPrt(String,' ',Work(LW4),nConf)
*         End If
*       End Do
        If ( NAC.eq.0 ) then
          ENER(1,ITER)=EMY
        Else
*         Do jRoot = 1,lRootSplit
*           ENER(jRoot,ITER) = Work(ipSplitE+jRoot-1)
            ENER(lRootSplit,ITER) = Work(ipSplitE+lRootSplit-1)
*           write(6,*)'ENER(jRoot,ITER)',ENER(jRoot,ITER)
            if (DBG) then
              write(6,*)'ITER :',ITER
              write(6,*)'ENER(lRootSplit,ITER)',ENER(lRootSplit,ITER)
            end if
*         End Do
        End if
************************************************************************
*            CLEANUP AFTER native DIAGONALIZATION                      *
************************************************************************
      call getmem('SplitV','FREE','Real',ipSplitV,iDimBlockA*iDimBlockA)
      call getmem('SplitE','FREE','Real',ipSplitE,iDimBlockA)
      CALL GETMEM('DHAM','FREE','REAL',ipDHAM,iDimBlockTri)
      CALL GETMEM('IREOTS','FREE','INTEGER',IREOTS,NAC)
      CALL GETMEM('AAblock','FREE','REAL',ipAABlock,nAAblock)
      CALL GETMEM('HONE','FREE','REAL',LOCONE,NAC**2)
      CALL GETMEM('IPCNF','FREE','INTE',LG1,NCNASM(LSYM))
      Call GetMem('iSel','Free','Integer',lSel,MxSpli)
      end if
*     ^ End of SplitCas/Native Hamiltonian diagonalization

      ELSE
*     ^... Do it IF there is CIRESTART
      iJOB=0
      Call f_Inquire('JOBOLD',Exist)
      If (Exist) iJOB=1
      If ( iJOB.eq.1) Write (6,*)
     &    ' Initial CI-vectors are read from JOBOLD'
      If ( iJOB.eq.0) Write (6,*)
     &    ' Initial CI-vectors are read from JOBIPH'
      If (iJOB.eq.1) Then
        If (JOBOLD.le.0) Then
          JOBOLD=20
          Call DaName(JOBOLD,'JOBOLD')
        End If
      Else
        JOBOLD=JOBIPH
      End If
      iDisk = 0
      Call IDafile(JOBOLD,2,iToc,15,iDisk)
      iDisk = iToc(4)
      Call GetMem('Scr1','Allo','Real',iTmp1,nConf)
*      Do i = 1,lRootSplit
        Call DDafile(JOBOLD,2,Work(iTmp1),nConf,iDisk)
        call GetMem('kcnf','allo','inte',ivkcnf,nactel)
        Call Reord2(NAC,NACTEL,LSYM,1,
     &              iWork(KICONF(1)),iWork(KCFTP),
     &              Work(iTmp1),C,iwork(ivkcnf))
        call GetMem('kcnf','free','inte',ivkcnf,nactel)
        Call Save_CI_vec(1,nConf,C,LuDavid)
*        Write (6,'(A,I2)') 'Start vector of root',i
*      if (DBG) then
        write(6,*)'LuDavid',LuDavid
        Write (String,'(A)') '(CI coefficient in CIRST)'
        Call dVcPrt(String,' ',C,nConf)
*      end if
*      End Do
      Call GetMem('Scr1','Free','Real',iTmp1,nConf)
      If (iJOB.eq.1) Then
        If(JOBOLD.gt.0.and.JOBOLD.NE.JOBIPH) Then
          Call DaClos(JOBOLD)
          JOBOLD=-1
        Else If(JOBOLD.gt.0) Then
          JOBOLD=-1
        End If
      End If

      END IF
*     ^End of do it IF there (is)/(isn't) CIRESTART

* -------------------------------------------------------------------- C
* -- CLEANUP AFTER SPLITCAS and DAVIDSON DIAGONALIZATION
* -------------------------------------------------------------------- C
*
*     LW4: TEMPORARY CI VECTOR IN CSF BASIS

      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
      IF (IPRLEV.GE.20) WRITE(6,1100) 'TERM_DAVID',LW4
      iDisk = IADR15(4)
*      Call Term_David(ICICH,ITERCI,lRootSplit,nConf,Work(LW4),JOBIPH,
*     &                LuDavid,iDisk)
*      write(6,*)'ITERCI :',ITERCI
      Call Term_David(ICICH,ITERCI,1,nConf,Work(LW4),JOBIPH,
     &                LuDavid,iDisk)
      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)

      if (DBG) then
        call cwtime(CSplitTot2,WSplitTot2)
        write(6,*) 'Total Time in SplitCAS:'
        write(6,*) 'CPU timing : ', CSplitTot2 - CSplitTot1
        write(6,*) 'W. timing  : ', WSplitTot2 - WSplitTot1
      end if
      Call qExit('SplitCTL')

      Return

1100  FORMAT(1X,/,1X,'WORK SPACE VARIABLES IN SUBR. CICTL: ',/,
     &       1X,'SUBSECTION: ',A,/,(1X,12I10,/))
      End
