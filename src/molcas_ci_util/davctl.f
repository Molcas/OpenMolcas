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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991, Jeppe Olsen                                      *
*               1991,1996, Markus P. Fuelscher                         *
*               2000, Thorstein Thorsteinsson                          *
************************************************************************
      Subroutine DavCtl(LW1, TUVX, IFINAL)
************************************************************************
*                                                                      *
*     CI control section                                               *
*                                                                      *
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
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated to use determinant based CI-procedures                 *
*       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
*     - updated for MOLCAS version 3                                   *
*       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*     - various modifications                                          *
*       T. Thorsteinsson, University of Lund, Sweden, 2000             *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension LW1(*), TUVX(*)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
#include "lucia_ini.fh"
      Parameter(ROUTINE='DAVCTL  ')
      Call qEnter('DAVCTL')
C
C -------------------------------------------------------------------- C
C -- INITIALIZE THE DAVIDSON DIAGONALIZATION
C -------------------------------------------------------------------- C
C
      lRoots=lRoots+hroots
      Call Ini_David(lRoots,nConf,nDet,nSel,nAc,LuDavid)

      IPRLEV=IPRLOC(3)
C
C -------------------------------------------------------------------- C
C --  COMPUTE THE DIAGONAL ELEMENTS OF THE HAMILTONIAN
C -------------------------------------------------------------------- C
C
C     LW4: TEMPORARY CI VECTOR IN CSF BASIS
C
      IF (IPRLEV.GE.20) WRITE(6,1100) 'INI_DAVID'
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
      IF (IPRLEV.GE.20) WRITE(6,1100) 'CIDIA',LW4
      IF (NAC. GT. 0)
     &       CALL CIDIA_CI_UTIL(NAC,NCONF,LSYM,WORK(LW4),LW1,TUVX,
     &                          LUDAVID)
C
C -------------------------------------------------------------------- C
C -- OBTAIN STARTING VECTORS
C -------------------------------------------------------------------- C
C
      mSel=nSel
      IF (NAC .NE. 0) THEN
         Call GetMem('iSel','Allo','Integer',lSel,mSel)
         Call GetMem('ExplE','Allo','Real',lExplE,mSel)
         Call GetMem('ExplV','Allo','Real',lExplV,mSel*mSel)
      ELSE
         lSel=1
         lExplE=1
         lExplV=1
      ENDIF
      IF (IPRLEV.GE.20 .AND. NAC .NE. 0)
     &        WRITE(6,1100) 'CSTART',LW4,lSel,lExplE,lExplV
      Call CStart_CI_Util(WORK(LW4),LW1,TUVX,
     &     iWork(lSel),Work(lExplE),Work(lExplV),IFINAL)

*MGD if nsel=nCSF_HEXS, save CI vectors
        If (N_ELIMINATED_GAS_MOLCAS.gt.0.and.
     &      nSel.ne.nConf.and.nSel.eq.nCSF_HEXS) then
          Do i = 1,lRoots
            Call Load_CI_vec(i,nConf,WORK(LW4),LuDavid)
            Call Save_tmp_CI_vec(i,nConf,WORK(LW4),LuDavid)
          End Do
        EndIf
      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)
C
C -------------------------------------------------------------------- C
C -- DIAGONALIZATION SECTION
C -------------------------------------------------------------------- C

* PAM Jun 2006: Gradual lowering of threshold in first 3 iterations
*      If ( Iter.gt.1 ) Threshold=THFACT*ABS(CONV(4,ITER-1))
C Energy threshold for CI convergence criterion:
      If ( Iter.eq.1 ) THEN
        Threshold=THREN
      ELSE IF (ITER.gt.1 .and. ITER.le.3) THEN
        ThrRule=THFACT*ABS(CONV(4,ITER-1))
        Threshold=(DBLE(4-ITER)*THREN + DBLE(ITER)*ThrRule)*0.25D0
      ELSE
        Threshold=THFACT*ABS(CONV(4,ITER-1))
      END IF
* End of new rule, PAM Jun 2006
      Threshold=Max(Threshold,1.0D-9)
      If(NAC.eq.0) then
        ESize=ABS(EMY)
      else
        ESize=ABS(Work(lExplE))
      end if
      Threshold=Max(Threshold,ESize*1.0D-14)

C     LW5: CONVERGENCE PARAMETERS
      Call GetMem('CI_conv','Allo','Real',LW5,2*lRoots*MAXJT)
      IF (IPRLEV.GE.20 .AND. NAC .NE. 0)
     &        WRITE(6,1100) 'DAVID',LW5,lSel,lExplE,lExplV
      ITERCI=1
      If ( NAC.eq.0 ) then
        ENER(1,ITER)=EMY
      Else
        If (( nSel.eq.nConf ).or.
     &     (N_ELIMINATED_GAS_MOLCAS.gt.0.and.nSel.eq.nCSF_HEXS)) then
          Do jRoot = 1,lRoots-hRoots
            ENER(jRoot,ITER) = Work(lExplE+jRoot-1)
          End Do
        Else
* PAM Jun 2006: Limit the number of CI iterations in the
* first few macroiterations:
          If(KTIGHT.eq.0) ItLimit=min(12*ITER,MAXJT)
          If(KTIGHT.eq.1) ItLimit=MAXJT
* PAM Oct 2006: Full precision if this is final CI.
          IF (ICIONLY.eq.1 .or. IFINAL.eq.2) THEN
            Threshold=Max(1.0D-9,ESize*1.0D-14)
            ITLIMIT=MAXJT
          END IF
* PAM Feb 2009: New code in david5.
*           Call David5(nAc,lSym,nDet,MAXJT,ITERCI,
*          Call David5(nAc,lSym,nDet,ItLimit,ITERCI,
*     &      Work(LW5),Threshold,LW1, TUVX,
*     &      iWork(lSel),Work(lExplE),Work(lExplV))
*
           Call David5(nDet,ItLimit,IterCI,Work(LW5),Threshold,
     &                 iWork(lSel),Work(lExplE),Work(lExplV),
     &                 LW1,TUVX)

          Do jRoot = 1,lRoots-hRoots
            ENER(jRoot,ITER) = Work(LW5+2*(jRoot-1)+
     &        2*(ITERCI-1)*lRoots)
          End Do
        End If
      End If
      Call GetMem('CI_conv','Free','Real',LW5,2*lRoots*MAXJT)
      IF (NAC .NE. 0) THEN
         Call GetMem('ExplV','Free','Real',lExplV,mSel*mSel)
         Call GetMem('ExplE','Free','Real',lExplE,mSel)
         Call GetMem('iSel','Free','Integer',lSel,mSel)
      ENDIF
      nSel = mSel
      lRoots=lRoots-hroots
*
C
C -------------------------------------------------------------------- C
C -- CLEANUP AFTER THE DAVIDSON DIAGONALIZATION
C -------------------------------------------------------------------- C
C
C     LW4: TEMPORARY CI VECTOR IN CSF BASIS
C
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
      IF (IPRLEV.GE.20) WRITE(6,1100) 'TERM_DAVID',LW4
      iDisk = IADR15(4)
      Call Term_David(ICICH,ITERCI,lRoots,nConf,Work(LW4),
     &                JOBIPH,LuDavid,iDisk)
      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)

      Call qExit('DAVCTL')

      Return

1100  FORMAT(1X,/,1X,'WORK SPACE VARIABLES IN SUBR. CICTL: ',/,
     &       1X,'SUBSECTION: ',A,/,(1X,12I10,/))
      End
      Subroutine Ini_David(nRoots,nConf,nDet,nSel,ntAsh,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Prepare address tables for further use by the Davidson           *
*     diagonalization scheme which is written in such a way that       *
*     the mechanism by which I/O is done is hidden in the subroutines  *
*     call load_xxx and save_xxx, where xxx stands may be one of       *
*     the following choices: H_diag, CI_vec, Sig_vec or Tmp_vec.       *
*     If possible (there is enough memory) a write through cache       *
*     mechanism is applied, that is to say all accessible memory is    *
*     used as a RAM-disk and dumped to physical disk in a FIFO mode.   *
*                                                                      *
*     calling arguments:                                               *
*     lRoots  : integer                                                *
*               number of roots to be optimized                        *
*     nConf   : integer                                                *
*               length of the CI vector in the CSF basis               *
*     nDet    : integer                                                *
*               length of the CI vector in the determinant basis       *
*     ntAsh   : integer                                                *
*               total number of active orbitals                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"

      Character*8 Label
      Real*8 Dum

      Call qEnter('Ini_David')

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Ini_David: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( nRoots.lt.0 ) then
         Write(6,*) 'Ini_David: nRoots less than zero'
         Write(6,*) 'nRoots = ',nRoots
         Call QTrace
         Call Abend
      Endif
      If ( nRoots.gt.mxRoot ) then
         Write(6,*) 'Ini_David: nRoots greater than mxRoot'
         Write(6,*) 'nRoots, mxRoot = ',nRoots, mxRoot
         Call QTrace
         Call Abend
      Endif
      If ( nDet.lt.0 ) then
         Write(6,*) 'Ini_David: nDet less than zero'
         Write(6,*) 'nDet = ',nDet
         Call QTrace
         Call Abend
      Endif
      If ( ntAsh.lt.0 ) then
         Write(6,*) 'Ini_David: ntAsh less than 0'
         Write(6,*) 'ntAsh = ',ntAsh
         Call QTrace
         Call Abend
      Endif
      If ( ntAsh.gt.mxAct ) then
         Write(6,*) 'Ini_David: ntAsh greater than mxAct'
         Write(6,*) 'ntAsh, mxAct = ',ntAsh, mxAct
         Call QTrace
         Call Abend
      Endif
      n_Roots=nRoots
*     Determine a reasonable nkeep
      nkeep=mxKeep*nRoots
      nkeep=min(nkeep,300)
      nkeep=max(nkeep,3*nRoots)
*
      istart=0
      nvec=nkeep

*     check the amount of available memory and decide which algorithm
*     is to be used to save intermediate results
      MxMemStk  = 0
      MxDiskStk = 0
      Call GetMem(' ','Max ','Real',Max_free_Mem,Max_free_Mem)
      Max_free_Mem = Max_free_Mem - 3*(nDet+4)
      Max_free_Mem = Max_free_Mem - 3*(nConf+4)
      Max_free_Mem = Max_free_Mem - 2*(ntAsh**3+4)
      Max_free_Mem = Max_free_Mem - 5*(ntAsh**2+4)
      Max_used_Mem = (1 + 2*nKeep+2*nRoots)*(nConf+4)
* Calculate how much memory is needed in the rest of the Davidson
      Memory_Needed = 0
      If (ntAsh .EQ. 0) Then
         Memory_Needed = 0
      Else If (nSel .EQ. nConf) Then
         Memory_Needed = 2*nSel + nSel*nSel
      Else
* First: davctl
         Memory_Needed = 2*nSel + nSel*nSel
* Now: david5
         lTmp1 = nKeep
         lTmp2 = lTmp1*lTmp1
         lTmp3 = (lTmp2+lTmp1)/2
         Memory_Needed = Memory_Needed + 5*nDet + lTmp1 +
     &                   3*lTmp2 + 2*lTmp3 + 3*nRoots*nSel
* Then: lucia_util
         Memory_Needed = Memory_Needed + Memory_Needed_Lucia
      End If
*
      If ( Max_free_Mem.lt.(nConf+4+Memory_Needed) ) then
        MxMemStk  = 0
        MxDiskStk = 1 + 2*nkeep + 2*nRoots
        save_mode = on_disk
      Else If ( Max_free_Mem.ge.Max_used_Mem + Memory_Needed ) then
        MxMemStk  = 1 + 2*nkeep + 2*nRoots
        MxDiskStk = 0
        save_mode = in_core
      Else
        MxMemStk  = Max_free_Mem/(nConf+4+Memory_Needed)
        MxDiskStk = 1 + 2*nkeep + 2*nRoots - mxMemStk
        save_mode = mixed_mode_2
        If ( mxMemStk.lt.(nkeep+1) ) save_mode = mixed_mode_1
      End If
CFUE  Call GetMem(' ','nFld',' ',nMemStk,nMemStk)
CFUE  nMemStk = nMemStk - 30
CFUE  If ( MxMemStk.gt.nMemStk ) then
CFUE    MxMemStk  = nMemStk
CFUE    MxDiskStk = 1 + 2*mxKeep*nRoots + 2*nRoots - mxMemStk
CFUE    save_mode = mixed_mode_2
CFUE    If ( mxMemStk.lt.(mxKeep*nRoots+1) ) save_mode = mixed_mode_1
CFUE  End If
      nMemStk = 0
      nDiskStk = 0

*     the diagonalization can be run in core:
*     allocate memory for all vectors that will be needed
      If ( save_mode.eq.in_core ) then
        H_diag_RecNo = RecNo((1),(1))
        Write(Label,'(A,I3.3)') 'HvRcN',H_diag_RecNo
        Call GetMem(Label,'Allo','Real',iMem,nConf)
        memory_address(H_diag_RecNo) = iMem
        CI_vec_RecNo = 0
        Do iRoot = 1,nkeep
          CI_vec_RecNo = RecNo((2),iRoot)
          Write(Label,'(A,I3.3)') 'CvRcN',CI_vec_RecNo
          Call GetMem(Label,'Allo','Real',iMem,nConf)
          memory_address(CI_vec_RecNo) = iMem
        End Do
        Sig_vec_RecNo = 0
        Do iRoot = 1,nKeep
          Sig_vec_RecNo = RecNo((3),iRoot)
          Write(Label,'(A,I3.3)') 'SvRcN',Sig_vec_RecNo
          Call GetMem(Label,'Allo','Real',iMem,nConf)
          memory_address(Sig_vec_RecNo) = iMem
        End Do
        Do iRoot = 1,nRoots
          tmp_CI_vec_RecNo = RecNo((4),iRoot)
          Write(Label,'(A,I3.3)') 'TmpCv',iRoot
          Call GetMem(Label,'Allo','Real',iMem,nConf)
          memory_address(tmp_CI_vec_RecNo) = iMem
        End Do
        Do iRoot = 1,nRoots
          tmp_Sig_vec_RecNo = RecNo((5),iRoot)
          Write(Label,'(A,I3.3)') 'TmpSv',iRoot
          Call GetMem(Label,'Allo','Real',iMem,nConf)
          memory_address(tmp_Sig_vec_RecNo) = iMem
        End Do
      End If

*     the diagonalization must be run out of core:
*     allocate disk space for all vectors that will be needed
      If ( save_mode.eq.on_disk ) then
        iDisk  = 0
        H_diag_RecNo = RecNo((1),(1))
        disk_address(H_diag_RecNo) = iDisk
        Call DDafile(LuDavid,0,Dum,nConf,iDisk)
        Do iRoot = 1,nkeep
          CI_vec_RecNo = RecNo((2),iRoot)
          disk_address(CI_vec_RecNo) = iDisk
          Call DDafile(LuDavid,0,Dum,nConf,iDisk)
        End Do
        Do iRoot = 1,nKeep
          Sig_vec_RecNo = RecNo((3),iRoot)
          disk_address(Sig_vec_RecNo) = iDisk
          Call DDaFile(LuDavid,0,Dum,nConf,iDisk)
        End Do
        Do iRoot = 1,nRoots
          tmp_CI_vec_RecNo = RecNo((4),iRoot)
          disk_address(tmp_CI_vec_RecNo) = iDisk
          Call DDaFile(LuDavid,0,Dum,nConf,iDisk)
        End Do
        Do iRoot = 1,nRoots
          tmp_Sig_vec_RecNo = RecNo((5),iRoot)
          disk_address(tmp_Sig_vec_RecNo) = iDisk
          Call DDaFile(LuDavid,0,Dum,nConf,iDisk)
        End Do
      End If

*     the diagonalization may be run in mixed mode:
*     allocate memory and disk space for all vectors that will be needed
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        Do nStk = 1,mxMemStk
          Write(Label,'(A,I4.4)') 'RAMD',nStk
          Call GetMem(Label,'Allo','Real',iMem,nConf)
          memory_address(nStk) = iMem
        End Do
        iDisk = 0
        Do nStk = 1,mxDiskStk
          disk_address(nStk) = iDisk
          Call DDaFile(LuDavid,0,Dum,nConf,iDisk)
        End Do
        Do nStk = 1,(mxMemStk+mxDiskStk)
          LblStk(nStk) = '                '
        End Do
        save_in_memory = .true.
      End If

      Call qExit('Ini_David')

      Return
      End
      Subroutine Term_David(ICICH,iter,lRoots,nConf,Vector,
     &                      JOBIPH,LuDavid,iDisk)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Terminate the Davidson diagonalization                           *
*                                                                      *
*     calling arguments:                                               *
*     ICICH   : integer                                                *
*               switch enabling root selection                         *
*     JOBIPH  : integer                                                *
*               logical unit number of the JOBIPH file                 *
*     iDisk   : integer                                                *
*               disk address of the first CI vector on JOBIPH          *
*     iter    : integer                                                *
*               iteration count of the final result                    *
*     nConf   : integer                                                *
*               length of the CI vector in the CSF basis               *
*     Vector  : array of real*8                                        *
*               temporary vector of length nConf                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Vector(nConf)


#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"

      Call qEnter('Term_David')

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Term_David: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iter.lt.0 ) then
         Write(6,*) 'Term_David: iter less than 0'
         Write(6,*) 'iter = ',iter
         Call QTrace
         Call Abend
      Endif
      If ( iter.gt.mxCiIt ) then
         Write(6,*) 'Term_David: iter greater than mxCiIt'
         Write(6,*) 'iter, mxCiIt = ',iter, mxCiIt
         Call QTrace
         Call Abend
      Endif

*     Restore the final CI vectors and save them for further use.
*     If the root selectioning option has been enabled calculate
*     also the overlap elemtents with the test vectors
      If ( ICICH.eq.1 ) then
        Call GetMem('CIovlp1','Allo','Real',lOvlp1,lRoots*lRoots)
        Call dCopy_(lRoots*lRoots,0.0d0, 0,Work(lOvlp1),(1))
        Call GetMem('CIovlp2','Allo','Real',lOvlp2,lRoots*lRoots)
        Call dCopy_(lRoots*lRoots,0.0d0, 0,Work(lOvlp2),(1))
      End If
      Do iRoot = 1,lRoots
        Call Load_tmp_CI_vec(iRoot,nConf,Vector,LuDavid)
        Call DDaFile(JOBIPH,1,Vector,nConf,iDisk)
        If ( ICICH.eq.1 ) then
          Call CIovlp(iRoot,Work(lOvlp1),Work(lOvlp2),Vector)
        End if
      End Do

*     If the root selectioning option has been enabled
*     make a new choice of the current roots
      If ( ICICH.eq.1 ) then
        Call CIselect(Work(lOvlp1),Work(lOvlp2))
        Call GetMem('CIovlp2','Free','Real',lOvlp2,lRoots*lRoots)
        Call GetMem('CIovlp1','Free','Real',lOvlp1,lRoots*lRoots)
      End If

*     deallocate memory which was used as records of the RAM disk
      If ( save_mode.ne.on_disk ) then
        Do iRecNo = 1,MxMemStk
          iMem = memory_address(iRecNo)
          Call GetMem(' ','Free','Real',iMem,nConf)
        End Do
      End If

      Call qExit('Term_David')

      Return
      End
      Subroutine Load_H_diag(nConf,H_diag,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load the diagonal approximation of the CI Hamiltonian for        *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     H_diag  : array of real*8                                        *
*               diagonal approximation of the CI Hamiltonian           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 H_diag(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Load_H_diag')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Load_H_diag: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy H_diag to new memory location
      If ( save_mode.eq.in_core ) then
        H_diag_RecNo = RecNo((1),(1))
        iMem = memory_address(H_diag_RecNo)
        Call dCopy_(nConf,Work(iMem),1,H_diag,1)
      End If

*     the diagonalization must be run out of core:
*     load H_diag from disk
      If ( save_mode.eq.on_disk ) then
        H_diag_RecNo = RecNo((1),(1))
        iDisk = disk_address(H_diag_RecNo)
        Call DDaFile(LuDavid,2,H_diag,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save H_diag
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A)') 'H_diag'
        Call page_in(KeyWord,nConf,H_diag,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Load_H_diag')

      Return
      End
      Subroutine Save_H_diag(nConf,H_diag,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save the diagonal approximation of the CI Hamiltonian for        *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     H_diag  : array of real*8                                        *
*               diagonal approximation of the CI Hamiltonian           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 H_diag(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Save_H_diag')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Save_H_diag: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy vector to new memory location
      If ( save_mode.eq.in_core ) then
        H_diag_RecNo = RecNo((1),(1))
        iMem = memory_address(H_diag_RecNo)
        Call dCopy_(nConf,H_diag,1,Work(iMem),1)
      End If

*     the diagonalization must be run out of core:
*     save H_diag on disk
      If ( save_mode.eq.on_disk ) then
        H_diag_RecNo = RecNo((1),(1))
        iDisk = disk_address(H_diag_RecNo)
        Call DDaFile(LuDavid,1,H_diag,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save H_diag
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A)') 'H_diag'
        Call page_out(KeyWord,nConf,H_diag,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Save_H_diag')

      Return
      End
      Subroutine Load_CI_vec(iRoot,nConf,CI_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load a CI vector                                                 *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     CI_vec  : array of real*8                                        *
*               CI vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 CI_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Load_CI_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Load_CI_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Load_CI_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.nkeep ) then
         Write(6,*) 'Load_CI_vec: iRoot greater than nkeep'
         Write(6,*) 'iRoot, nkeep = ',iRoot, nkeep
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the CI vector to new memory location
      If ( save_mode.eq.in_core ) then
        CI_vec_RecNo = RecNo((2),iRoot)
        iMem = memory_address(CI_vec_RecNo)
        Call dCopy_(nConf,Work(iMem),1,CI_vec,1)
      End If

*     the diagonalization must be run out of core:
*     load the CI vector from disk
      If ( save_mode.eq.on_disk ) then
        CI_vec_RecNo = RecNo((2),iRoot)
        iDisk = disk_address(CI_vec_RecNo)
        Call DDaFile(LuDavid,2,CI_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the CI vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        CI_vec_PageNo = PageNo(iRoot)
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'CI_vec',CI_vec_PageNo
        Call page_in(KeyWord,nConf,CI_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Load_CI_vec')

      Return
      End
      Subroutine Save_CI_vec(iRoot,nConf,CI_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save a CI vector                                                 *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     CI_vec  : array of real*8                                        *
*               CI vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 CI_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Save_CI_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Save_CI_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Save_CI_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.nkeep ) then
         Write(6,*) 'Save_CI_vec: iRoot greater than nkeep'
         Write(6,*) 'iRoot, nkeep = ',iRoot, nkeep
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the CI vector to new memory location
      If ( save_mode.eq.in_core ) then
        CI_vec_RecNo = RecNo((2),iRoot)
        iMem = memory_address(CI_vec_RecNo)
        Call dCopy_(nConf,CI_vec,1,Work(iMem),1)
      End If

*     the diagonalization must be run out of core:
*     save the CI vector on disk
      If ( save_mode.eq.on_disk ) then
        CI_vec_RecNo = RecNo((2),iRoot)
        iDisk = disk_address(CI_vec_RecNo)
        Call DDaFile(LuDavid,1,CI_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the CI vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        CI_vec_PageNo = PageNo(iRoot)
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'CI_vec',CI_vec_PageNo
        Call page_out(KeyWord,nConf,CI_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Save_CI_vec')

      Return
      End
      Subroutine Load_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load a sigma vector                                              *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Sig_vec : array of real*8                                        *
*               sigma vector                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Sig_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Load_Sig_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Load_Sig_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Load_Sig_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.nkeep ) then
         Write(6,*) 'Load_Sig_vec: iRoot greater than nkeep'
         Write(6,*) 'iRoot, nkeep = ',iRoot, nkeep
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the sigma vector to new memory location
      If ( save_mode.eq.in_core ) then
        Sig_vec_RecNo = RecNo((3),iRoot)
        iMem = memory_address(Sig_vec_RecNo)
        Call dCopy_(nConf,Work(iMem),1,Sig_vec,1)
      End If

*     the diagonalization must be run out of core:
*     load the sigma vector from disk
      If ( save_mode.eq.on_disk ) then
        Sig_vec_RecNo = RecNo((3),iRoot)
        iDisk = disk_address(Sig_vec_RecNo)
        Call DDaFile(LuDavid,2,Sig_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the sigma vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        Sig_vec_PageNo = PageNo(iRoot)
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'Sig_vec',Sig_vec_PageNo
        Call page_in(KeyWord,nConf,Sig_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Load_Sig_vec')

      Return
      End
      Subroutine Save_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save a sigma vector                                              *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Sig_vec : array of real*8                                        *
*               sigma vector                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Sig_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Save_Sig_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Save_Sig_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Save_Sig_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.nkeep ) then
         Write(6,*) 'Save_Sig_vec: iRoot greater than nkeep'
         Write(6,*) 'iRoot, nkeep = ',iRoot, nkeep
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the sigma vector to new memory location
      If ( save_mode.eq.in_core ) then
        Sig_vec_RecNo = RecNo((3),iRoot)
        iMem = memory_address(Sig_vec_RecNo)
        Call dCopy_(nConf,Sig_vec,1,Work(iMem),1)
      End If

*     the diagonalization must be run out of core:
*     save the sigma vector on disk
      If ( save_mode.eq.on_disk ) then
        Sig_vec_RecNo = RecNo((3),iRoot)
        iDisk = disk_address(Sig_vec_RecNo)
        Call DDaFile(LuDavid,1,Sig_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the sigma vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        Sig_vec_PageNo = PageNo(iRoot)
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'Sig_vec',Sig_vec_PageNo
        Call page_out(KeyWord,nConf,Sig_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Save_Sig_vec')

      Return
      End
      Subroutine Load_tmp_CI_vec(iRoot,nConf,CI_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load a temporary CI vector                                       *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     CI_vec  : array of real*8                                        *
*               CI vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 CI_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Load_tmp_CI_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Load_tmp_CI_vec: nConf less than'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Load_tmp_CI_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.n_Roots ) then
         Write(6,*) 'Load_tmp_CI_vec: iRoot greater than nRoots'
         Write(6,*) 'iRoot, nRoots = ',iRoot, n_Roots
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the CI vector to new memory location
      If ( save_mode.eq.in_core ) then
        tmp_CI_vec_RecNo = RecNo((4),iRoot)
        iMem = memory_address(tmp_CI_vec_RecNo)
        Call dCopy_(nConf,Work(iMem),1,CI_vec,1)
      End If

*     the diagonalization must be run out of core:
*     load the CI vector from disk
      If ( save_mode.eq.on_disk ) then
        tmp_CI_vec_RecNo = RecNo((4),iRoot)
        iDisk = disk_address(tmp_CI_vec_RecNo)
        Call DDaFile(LuDavid,2,CI_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the CI vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'tmp_CI_vec',iRoot
        Call page_in(KeyWord,nConf,CI_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Load_tmp_CI_vec')

      Return
      End
      Subroutine Save_tmp_CI_vec(iRoot,nConf,CI_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save a temporary CI vector                                       *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     CI_vec  : array of real*8                                        *
*               CI vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 CI_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Save_tmp_CI_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Save_tmp_CI_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Save_tmp_CI_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.n_Roots ) then
         Write(6,*) 'Save_tmp_CI_vec: iRoot greater than nRoots'
         Write(6,*) 'iRoot, nRoots = ',iRoot, n_Roots
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the CI vector to new memory location
      If ( save_mode.eq.in_core ) then
        tmp_CI_vec_RecNo = RecNo((4),iRoot)
        iMem = memory_address(tmp_CI_vec_RecNo)
        Call dCopy_(nConf,CI_vec,1,Work(iMem),1)
      End If

*     the diagonalization must be run out of core:
*     save the CI vector on disk
      If ( save_mode.eq.on_disk ) then
        tmp_CI_vec_RecNo = RecNo((4),iRoot)
        iDisk = disk_address(tmp_CI_vec_RecNo)
        Call DDaFile(LuDavid,1,CI_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the CI vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'tmp_CI_vec',iRoot
        Call page_out(KeyWord,nConf,CI_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Save_tmp_CI_vec')

      Return
      End
      Subroutine Load_tmp_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load a temporary sigma vector                                    *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Sig_vec : array of real*8                                        *
*               sigma vector                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Sig_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Load_tmp_Sig_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Load_tmp_Sig_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Load_tmp_Sig_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.n_Roots ) then
         Write(6,*) 'Load_tmp_Sig_vec: iRoot greater than nRoots'
         Write(6,*) 'iRoot = ',n_Roots
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the sigma vector to new memory location
      If ( save_mode.eq.in_core ) then
        tmp_Sig_vec_RecNo = RecNo((5),iRoot)
        iMem = memory_address(tmp_Sig_vec_RecNo)
        Call dCopy_(nConf,Work(iMem),1,Sig_vec,1)
      End If

*     the diagonalization must be run out of core:
*     load the sigma vector from disk
      If ( save_mode.eq.on_disk ) then
        tmp_Sig_vec_RecNo = RecNo((5),iRoot)
        iDisk = disk_address(tmp_Sig_vec_RecNo)
        Call DDaFile(LuDavid,2,Sig_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the sigma vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'tmp_Sig_vec',iRoot
        Call page_in(KeyWord,nConf,Sig_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Load_tmp_Sig_vec')

      Return
      End
      Subroutine Save_tmp_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save a temporary sigma vector                                    *
*     further use by the Davidson diagonalization scheme               *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Sig_vec : array of real*8                                        *
*               sigma vector                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Sig_vec(nConf)


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"

      Character*16 KeyWord

      Call qEnter('Save_tmp_Sig_vec')
      Call Timing(WTC_1,Swatch,Swatch,Swatch)

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'Save_tmp_Sig_vec: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.lt.0 ) then
         Write(6,*) 'Save_tmp_Sig_vec: iRoot less than 0'
         Write(6,*) 'iRoot = ',iRoot
         Call QTrace
         Call Abend
      Endif
      If ( iRoot.gt.n_Roots ) then
         Write(6,*) 'Save_tmp_Sig_vec: iRoot greater than nRoots'
         Write(6,*) 'iRoot, nRoots = ',iRoot, n_Roots
         Call QTrace
         Call Abend
      Endif

*     the diagonalization can be run in core:
*     copy the sigma vector to new memory location
      If ( save_mode.eq.in_core ) then
        tmp_Sig_vec_RecNo = RecNo((5),iRoot)
        iMem = memory_address(tmp_Sig_vec_RecNo)
        Call dCopy_(nConf,Sig_vec,1,Work(iMem),1)
      End If

*     the diagonalization must be run out of core:
*     save the sigma vector on disk
      If ( save_mode.eq.on_disk ) then
        tmp_Sig_vec_RecNo = RecNo((5),iRoot)
        iDisk = disk_address(tmp_Sig_vec_RecNo)
        Call DDaFile(LuDavid,1,Sig_vec,nConf,iDisk)
      End If

*     the diagonalization may be run in mixed mode:
*     use the write through cache mechanism to load and save
*     the sigma vector
      If ( save_mode.eq.mixed_mode_1 .or.
     &     save_mode.eq.mixed_mode_2      ) then
        KeyWord = '                '
        Write(KeyWord,'(A,I4.4)') 'tmp_Sig_vec',iRoot
        Call page_out(KeyWord,nConf,Sig_vec,LuDavid)
      End If

      Call Timing(WTC_2,Swatch,Swatch,Swatch)
      WTC_2 = WTC_2 - WTC_1
      WTC_3 = WTC_3 + WTC_2
      Call qExit('Save_tmp_Sig_vec')

      Return
      End
      Subroutine page_out(KeyWord,nConf,Vector,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Save any vector for further use by the Davidson diagonalization  *
*     Labels identifying the vectors are kept in a stack and to        *
*     minimize a write through cache strategy is applied               *
*                                                                      *
*     calling arguments:                                               *
*     KeyWord : character*16                                           *
*               record identifier                                      *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Vector  : array of real*8                                        *
*               any vector of length nConf                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Vector(nConf)
      Character*16 KeyWord


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"

      Call qEnter('page_out')

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'page_out: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif

*     serach for a matching record identifier
      nStk = 0
      Do iStk = 1,(mxMemStk+mxDiskStk)
        If ( LblStk(iStk).eq.KeyWord ) nStk = iStk
      End Do

*     there is a matching record identifier:
*     overwrite the current record
      If ( nStk.ne.0 ) then
        If ( nStk.le.mxMemStk ) then
          iMem = memory_address(nStk)
          Call dCopy_(nConf,Vector,1,Work(iMem),1)
        Else
          iDisk = disk_address(nStk-mxMemStk)
          Call DDaFile(LuDavid,1,Vector,nConf,iDisk)
        End If
      End if

*     there is no matching record identifier:
*     create a new record
      If ( nStk.eq.0 ) then
        If ( save_mode.eq.mixed_mode_1 ) then
          If ( KeyWord(1:6).eq.'CI_vec' ) then
            If ( save_in_memory ) then
              nMemStk = nMemStk + 1
              iMem = memory_address(nMemStk)
              Call dCopy_(nConf,Vector,1,Work(iMem),1)
              LblStk(nMemStk) = KeyWord
              if ( nMemStk.eq.mxMemStk ) save_in_memory = .false.
            Else
              nMemStk = nMemStk + 1
              If ( nMemStk.gt.mxMemStk ) nMemStk = 1
              iMem = memory_address(nMemStk)
              nDiskStk = nDiskStk + 1
              If ( nDiskStk.gt.mxDiskStk ) nDiskStk = 1
              iDisk = disk_address(nDiskStk)
              Call DDaFile(LuDavid,1,Work(iMem),nConf,iDisk)
              Call dCopy_(nConf,Vector,1,Work(iMem),1)
              LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
              LblStk(nMemStk) = KeyWord
            End If
          Else
            nDiskStk = nDiskStk + 1
            If ( nDiskStk.gt.mxDiskStk ) nDiskStk = 1
            iDisk = disk_address(nDiskStk)
            Call DDaFile(LuDavid,1,Vector,nConf,iDisk)
            LblStk(mxMemStk+nDiskStk) = KeyWord
          End If
        End If
        If ( save_mode.eq.mixed_mode_2 ) then
          If ( save_in_memory ) then
            nMemStk = nMemStk + 1
            iMem = memory_address(nMemStk)
            Call dCopy_(nConf,Vector,1,Work(iMem),1)
            LblStk(nMemStk) = KeyWord
            If ( nMemStk.eq.mxMemStk ) save_in_memory = .false.
          Else
            nMemStk = nMemStk + 1
            If ( nMemStk.gt.mxMemStk ) nMemStk = 1
            iMem = memory_address(nMemStk)
            nDiskStk = nDiskStk + 1
            If ( nDiskStk.gt.mxDiskStk ) nDiskStk = 1
            iDisk = disk_address(nDiskStk)
            Call DDaFile(LuDavid,1,Work(iMem),nConf,iDisk)
            Call dCopy_(nConf,Vector,1,Work(iMem),1)
            LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
            LblStk(nMemStk) = KeyWord
          End if
        End If
      End if

      Call qExit('page_out')

      Return
      End
      Subroutine page_in(KeyWord,nConf,Vector,LuDavid)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Load any vector for further use by the Davidson diagonalization  *
*     which has been saved by the write through cache mechanism        *
*                                                                      *
*     calling arguments:                                               *
*     KeyWord : character*16                                           *
*               record identifier                                      *
*     nConf   : integer                                                *
*               length of the vector H_diag                            *
*     Vector  : array of real*8                                        *
*               any vector of length nConf                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

      Real*8 Vector(nConf)
      Character*16 KeyWord


#include "rasdim.fh"

#include "davctl.fh"
#include "WrkSpc.fh"

      Call qEnter('page_in')

*     check input arguments
      If ( nConf.lt.0 ) then
         Write(6,*) 'page_in: nConf less than 0'
         Write(6,*) 'nConf = ',nConf
         Call QTrace
         Call Abend
      Endif

*     serach for a metching record identifier
      nStk = 0
      Do iStk = 1,(mxMemStk+mxDiskStk)
        If ( LblStk(iStk).eq.KeyWord ) nStk = iStk
      End Do
      If ( nStk.eq.0 ) then
         Write(6,*) 'page_in: nStk equal 0'
         Write(6,*) 'nStk = ',nStk
         Call QTrace
         Call Abend
      Endif

      If ( nStk.le.mxMemStk ) then
        iMem = memory_address(nStk)
        Call dCopy_(nConf,Work(iMem),1,Vector,1)
      Else
        iDisk = disk_address(nStk-mxMemStk)
        Call DDaFile(LuDavid,2,Vector,nConf,iDisk)
      End If

      Call qExit('page_in')

      Return
      End
      Integer Function PageNo(iRoot)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Compute the page number of a vector                              *
*                                                                      *
*     calling arguments:                                               *
*     iRoot   : integer                                                *
*               root number                                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)


#include "rasdim.fh"

#include "davctl.fh"

*     Call qEnter('PageNo')

      itmp1 = iRoot
      If (iRoot.gt.n_Roots) then
        itmp1=n_Roots+mod(istart+iRoot-n_Roots-1,nvec-n_Roots)+1
      EndIf
      PageNo = itmp1

*     Call qExit('PageNo')

      Return
      End
      Integer Function RecNo(itype,iRoot)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Compute the Record number of a vector                            *
*                                                                      *
*     calling arguments:                                               *
*     itype   : integer                                                *
*               vector type: 1 = H_diag                                *
*                            2 = CI_vec                                *
*                            3 = Sig_vec                               *
*     iRoot   : integer                                                *
*               root number                                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)


#include "rasdim.fh"

#include "davctl.fh"

*     Call qEnter('RecNo')

      RecNo = 0
      If ( itype.eq.1 ) then
        H_diag_RecNo = 1
        RecNo = H_diag_RecNo
      Else If ( itype.eq.2 ) then
        CI_vec_RecNo = 1+PageNo(iRoot)
        RecNo = CI_vec_RecNo
      Else If ( itype.eq.3 ) then
        Sig_vec_RecNo = 1+nkeep+PageNo(iRoot)
        RecNo = Sig_vec_RecNo
      Else If ( itype.eq.4 ) then
        tmp_CI_vec_RecNo = 1+2*nKeep+iRoot
        RecNo = tmp_CI_vec_RecNo
      Else If ( itype.eq.5 ) then
        tmp_Sig_vec_RecNo = 1+2*nKeep+n_Roots+iRoot
        RecNo = tmp_Sig_vec_RecNo
      Else
        Write(6,*) 'RecNo: itype does not match'
        Write(6,*) 'itype = ',itype
        Call QTrace
        Call Abend
      End If

*     Call qExit('RecNo')

      Return
      End
