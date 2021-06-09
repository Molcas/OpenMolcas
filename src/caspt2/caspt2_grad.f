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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine GrdIni
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
#include "pt2_guga.fh"
C
      !! Allocate lagrangian terms
      nBasTr = 0
      nBasSq = 0
      nOLag = 0
      nCLag = 0
      DO iSym = 1, nSym
        nBasI = nBas(iSym)
        nBasTr = nBasTr + nBasI*(nBasI+1)/2
        nBasSq = nBasSq + nBasI*nBasI
        nCLag = nCLag + nState*nCSF(iSym)
      END DO
      nOLag = nBasSq
      nSLag = nState*(nState-1)/2
      nWLag = nBasSq
C
      Call GETMEM('DPT2   ','ALLO','REAL',ipDPT2   ,nBasSq)
      Call GETMEM('DPT2C  ','ALLO','REAL',ipDPT2C  ,nBasSq)
      Call GETMEM('DPT2AO ','ALLO','REAL',ipDPT2AO ,nBasTr)
      Call GETMEM('DPT2CAO','ALLO','REAL',ipDPT2CAO,nBasTr)
C
      CALL GETMEM('CLAG   ','ALLO','REAL',ipCLag   ,nCLag)
      CALL GETMEM('OLAG   ','ALLO','REAL',ipOLag   ,nOLag)
      CALL GETMEM('SLAG   ','ALLO','REAL',ipSLag   ,nSLag)
      CALL GETMEM('WLAG   ','ALLO','REAL',ipWLag   ,nWLag)
C     write(6,*) "nclag,nolag,nslag"
C     write(6,*)  nclag, nolag, nslag
C     write(6,*) ipclag,ipolag,ipslag
      Call DCopy_(nBasSq,0.0D+00,0,Work(ipDPT2)   ,1)
      Call DCopy_(nBasSq,0.0D+00,0,Work(ipDPT2C)  ,1)
      Call DCopy_(nBasTr,0.0D+00,0,Work(ipDPT2AO) ,1)
      Call DCopy_(nBasTr,0.0D+00,0,Work(ipDPT2CAO),1)
C
      Call DCopy_(nCLag ,0.0D+00,0,Work(ipCLag),1)
      Call DCopy_(nOLag ,0.0D+00,0,Work(ipOLag),1)
      Call DCopy_(nSLag ,0.0D+00,0,Work(ipSLag),1)
      Call DCopy_(nWLag ,0.0D+00,0,Work(ipWLag),1)
C
      !! LuGamma should be 60, but this record is used in MCLR, so
      !! have to use a different value. This number has to be consistent
      !! with that in ALASKA (integral_util/prepp.f).
C     LuGamma = 60
C
      Return
C
      End Subroutine GrdIni
C
C-----------------------------------------------------------------------
C
      Subroutine GrdCls
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Character(Len=16) mstate1
      LOGICAL DEBUG
C
      DEBUG = .FALSE.
      Call Molcas_Open(LuPT2,'PT2_Lag')
C     Write (LuPT2,*) BSHIFT
      !! configuration Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "CLag"
      Do i = 1, nCLag
C       if (abs(work(ipclag+i-1)).le.1.0d-10) work(ipclag+i-1)=0.0d+00
        Write (LuPT2,*) Work(ipClag+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipClag+i-1)
      End Do
      !! orbital Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "OLag"
      Do i = 1, nOLag
C       if (abs(work(ipolag+i-1)).le.1.0d-10) work(ipolag+i-1)=0.0d+00
        Write (LuPT2,*) Work(ipOlag+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipOlag+i-1)
      End Do
      !! state Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "SLag"
      Do i = 1, nSLag
        Write (LuPT2,*) Work(ipSlag+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipSlag+i-1)
      End Do
C
C
C
      !! renormalization contributions (read in OUT_PT2)
      If (DEBUG) write(6,*) "WLag"
      Do i = 1, nbast*(nbast+1)/2 !! nWLag
        Write (LuPT2,*) Work(ipWlag+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipWlag+i-1)
      End Do
C     write(6,*) "dpt2"
      !! D^PT2 in MO (read in OUT_PT2)
      If (DEBUG) write(6,*) "DPT2"
      Do i = 1, nBasSq
C       write(6,*) i,work(ipdpt2+i-1)
        Write (LuPT2,*) Work(ipDPT2+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipDPT2+i-1)
      End Do
C     write(6,*) "dpt2c"
      !! D^PT2(C) in MO (read in OUT_PT2)
      If (DEBUG) write(6,*) "DPT2C"
      Do i = 1, nBasSq
C       write(6,*) i,work(ipdpt2c+i-1)
        Write (LuPT2,*) Work(ipDPT2C+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipDPT2C+i-1)
      End Do
C
C
C
      !! D^PT2 in AO (not used?)
      Do i = 1, nBasTr
        Write (LuPT2,*) Work(ipDPT2AO+i-1)
      End Do
      !! D^PT2(C) in AO (not used?)
      Do i = 1, nBasTr
        Write (LuPT2,*) Work(ipDPT2CAO+i-1)
      End Do
      Close (LuPT2)
C
      Call GETMEM('DPT2   ','FREE','REAL',ipDPT2   ,nBasSq)
      Call GETMEM('DPT2C  ','FREE','REAL',ipDPT2C  ,nBasSq)
      Call GETMEM('DPT2AO ','FREE','REAL',ipDPT2AO ,nBasTr)
      Call GETMEM('DPT2CAO','FREE','REAL',ipDPT2CAO,nBasTr)
C
      CALL GETMEM('CLAG   ','FREE','REAL',ipCLag   ,nCLag)
      CALL GETMEM('OLAG   ','FREE','REAL',ipOLag   ,nOLag)
      CALL GETMEM('SLAG   ','FREE','REAL',ipSLag   ,nSLag)
      CALL GETMEM('WLAG   ','FREE','REAL',ipWLag   ,nWLag)
C
      !! Prepare for MCLR
C     If (Method.eq.'CASPT2  ') Then
         iGo = 0
         Call Put_iScalar('SA ready',iGo)
         mstate1 = '****************'
         Call Put_cArray('MCLR Root',mstate1,16)
C     End If
C       write(6,*) "5"
C     write(6,*) "LuGamma is ", LuGamma
C     write(6,*) "bshift =", bshift
C     Call Put_dScalar('BSHIFT',BSHIFT)
C
C
      Return
C
      End Subroutine GrdCls
C
C-----------------------------------------------------------------------
C
      Subroutine ModDip
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      CALL GETMEM('DMs1   ','ALLO','REAL',ipDMs1,3*nRoots)
      CALL GETMEM('DMs2   ','ALLO','REAL',ipDMs2,3*lRoots)
      Call Get_dArray('Last Dipole Moments',Work(ipDMs2),3*LROOTS)
      Do i = 1, lRoots
        j = Root2State(i)
        If (j.eq.0) Cycle
        Call DCopy_(3,Work(ipDMs2+3*(i-1)),1,Work(ipDMs1+3*(j-1)),1)
      End Do
      Call Put_dArray('Last Dipole Moments',Work(ipDMs1),3*nROOTS)
      CALL GETMEM('DMs1   ','FREE','REAL',ipDMs1,3*nRoots)
      CALL GETMEM('DMs2   ','FREE','REAL',ipDMs2,3*lRoots)
C
      Return
C
      End Subroutine ModDip
