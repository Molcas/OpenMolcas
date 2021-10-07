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
      !! nStLag is the number of states involved in the Lagrangian
      !! nStLag = nState for (X)MS/XDW/RMS-CASPT2
      !!        = 1      for SS-CASPT2
      If (IFMSCOUP) Then
        nStLag    = nState
      Else
        nStLag    = 1
      End If
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
      nSLag = nState*nState
      nWLag = nBasSq
C
      Call GETMEM('DPT2   ','ALLO','REAL',ipDPT2   ,nBasSq)
      Call GETMEM('DPT2C  ','ALLO','REAL',ipDPT2C  ,nBasSq)
      Call GETMEM('DPT2AO ','ALLO','REAL',ipDPT2AO ,nBasTr)
      Call GETMEM('DPT2CAO','ALLO','REAL',ipDPT2CAO,nBasTr)
C
      CALL GETMEM('CLAG   ','ALLO','REAL',ipCLag   ,nCLag*2)
      CALL GETMEM('OLAG   ','ALLO','REAL',ipOLag   ,nOLag*2)
      CALL GETMEM('SLAG   ','ALLO','REAL',ipSLag   ,nSLag)
      CALL GETMEM('WLAG   ','ALLO','REAL',ipWLag   ,nWLag)
C     write(6,*) "nclag,nolag,nslag"
C     write(6,*)  nclag, nolag, nslag
C     write(6,*) ipclag,ipolag,ipslag
      Call DCopy_(nBasSq ,[0.0D+00],0,Work(ipDPT2)   ,1)
      Call DCopy_(nBasSq ,[0.0D+00],0,Work(ipDPT2C)  ,1)
      Call DCopy_(nBasTr ,[0.0D+00],0,Work(ipDPT2AO) ,1)
      Call DCopy_(nBasTr ,[0.0D+00],0,Work(ipDPT2CAO),1)
C
      Call DCopy_(nCLag*2,[0.0D+00],0,Work(ipCLag),1)
      Call DCopy_(nOLag*2,[0.0D+00],0,Work(ipOLag),1)
      Call DCopy_(nSLag  ,[0.0D+00],0,Work(ipSLag),1)
      Call DCopy_(nWLag  ,[0.0D+00],0,Work(ipWLag),1)
C
      CALL GETMEM('FIFA   ','ALLO','REAL',ipFIFA   ,nBasSq)
      CALL GETMEM('FIMO   ','ALLO','REAL',ipFIMO   ,nBasSq)
C
      !! FIFASA is constructed with state-averaged density always
      !! FIFA   can be state-specific or dynamically weighted
      !! FIMO   is uniquely determined, but the basis can be
      !!        either natural or quasi-canonical
      If (IFXMS) Then
        CALL GETMEM('FIFASA ','ALLO','REAL',ipFIFASA ,nBasSq)
        norbi=norb(1)
      End If
C
      If (IFDW) Then
        CALL GETMEM('OMGDER ','ALLO','REAL',ipOMGDER ,nState**2)
        Call DCopy_(nState**2,[0.0D+00],0,Work(ipOMGDER),1)
      End If
C
      !! LuGamma should be 60, but this record is used in MCLR, so
      !! have to use a different value. This number has to be consistent
      !! with that in ALASKA (integral_util/prepp.f).
C     LuGamma = 60
C
      !! Some Lagrangians for each state are constructed in ipCLag or
      !! ipOLag. The full (sum over all states, in particular for
      !! MS-CASPT2) configuration and orbital Lagrangians are then
      !! constructed in ipOLagFull. For SS-CASPT2, ipCLag and
      !! ipCLagFull, for instance, will be identical.
      ipCLagFull = ipCLag + nCLag
      ipOLagFull = ipOLag + nOLag
C
      !! Whether the Fock matrix (eigenvalues) is constructed with
      !! the state-averaged density matrix or not.
      !! The name of the variable is like state-specific DM,
      !! but not necessarily state-specific. It is a matter of the
      !! structure of WORK(LDWGT) array or matrix.
      !! WORK(LDWGT) is a matrix form for SS- and MS-CASPT2 with
      !! state-specific DM, XDW-CASPT2, and RMS-CASPT2, while it is an
      !! array for SS- and MS-CASPT2 with state-averaged DM (with SADREF
      !! option) and XMS-CASPT2.
      If (IFSADREF.or.nState.eq.1.or.(IFXMS.and..not.IFDW)) Then
        IFSSDM = .false.
      Else
        IFSSDM = .true.
      End If
C
      If (isNAC) Then
        If (isCSF) Then
          Write (6,'(" Compute derivative coupling vectors for ",
     *               2i3)') iRoot1,iRoot2
        Else
          Write (6,'(" Compute interstate coupling vectors for ",
     *               2i3)') iRoot1,iRoot2
        End If
        Write (6,*)
        If (isCSF) Then
          Write (6,'(" Compute the full derivative coupling",
     *               " (include the CSF contribution)")')
          Write (6,'(" Note that CI and CSF derivative couplings are",
     *               " not strictly separable,")')
          Write (6,'(" because the CSF derivative coupling is",
     *               " dependent on wavefunction parameters,")')
          Write (6,'(" are the response is included in the CI",
     *               " derivative coupling.")')
          Write (6,'(" In addition, derivative coupling vectors are",
     *               " not translationally invariant.")')
        Else
          Write (6,'(" The CSF derivative coupling is not computed")')
        End If
        Write (6,*)
C
        Call GETMEM('DPT2Canti','ALLO','REAL',ipDPT2Canti,nBasSq)
        Call DCopy_(nBasSq ,[0.0D+00],0,Work(ipDPT2Canti),1)
      End If
C
      Return
C
      End Subroutine GrdIni
C
C-----------------------------------------------------------------------
C
      Subroutine GrdCls(UEFF,U0,H0)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Dimension UEFF(nState,nState),U0(nState,nState),H0(nState,nState)
      Character(Len=16) mstate1
      LOGICAL DEBUG
C
      Dimension HEFF1(nState,nState),WRK1(nState,nState),
     *          WRK2(nState,nState)
C
      !! Add XMS specific terms
      !! Note that ipCLagFull is in natural CSF basis,
      !! so everything in this subroutine has to be done in natural
      If (IFDW) Then
        !! Construct Heff[1] in XMS basis
        Call DCopy_(nState*nState,[0.0D+00],0,HEFF1,1)
        Do ilStat = 1, nState
         HEFF1(ilStat,ilStat) = REFENE(ilStat)
        End Do
        Call DGEMM_('T','N',nState,nState,nState,
     *              1.0D+00,U0,nState,HEFF1,nState,
     *              0.0D+00,WRK1,nState)
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,WRK1,nState,U0,nState,
     *              0.0D+00,HEFF1,nState)
C
        !! Derivative of Heff[1] in XMS basis
        !! It is transformed with U0, so the contribution has to be
        !! considered when we construct the auxiliary density in the
        !! XMS-specific term
        Do ilStat = 1, nState
          iloc = ilStat+nState*(ilStat-1)
          Ebeta = HEFF1(ilStat,ilStat)
          Do jlStat = 1, nState
            jloc = jlStat+nState*(jlStat-1)
            Ealpha = HEFF1(jlStat,jlStat)
            OMGDER = Work(ipOMGDER+ilStat-1+nState*(jlStat-1))
            Factor = 0.0d+00
            Do klStat = 1, nState
              Egamma = HEFF1(klStat,klStat)
              Factor = Factor + exp(-zeta*(Ealpha - Egamma)**2)
            End Do
            !! derivative of alpha-beta
            DERAB = -ZETA*EXP(-ZETA*(Ealpha-Ebeta)**2)/Factor
            Scal = 2.0D+00*DERAB*(Ealpha-Ebeta)*OMGDER
            Work(ipSLag+jloc-1) = Work(ipSLag+jloc-1) + Scal
            Work(ipSLag+iloc-1) = Work(ipSLag+iloc-1) - Scal
            !! derivative of alpha-gamma
            Do klStat = 1, nState
              kloc = klStat+nState*(klStat-1)
              Egamma = HEFF1(klStat,klStat)
              DERAC = ZETA*EXP(-ZETA*(Ealpha-Ebeta)**2)/(Factor*Factor)
     *                    *EXP(-ZETA*(Ealpha-Egamma)**2)
              Scal = 2.0D+00*DERAC*(Ealpha-Egamma)*OMGDER
              Work(ipSLag+jloc-1) = Work(ipSLag+jloc-1) + Scal
              Work(ipSLag+kloc-1) = Work(ipSLag+kloc-1) - Scal
            End Do
          End Do
        End Do
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,U0,nState,Work(ipSLag),nState,
     *              0.0D+00,WRK2,nState)
        Call DGEMM_('N','T',nState,nState,nState,
     *              1.0D+00,WRK2,nState,U0,nState,
     *              0.0D+00,WRK1,nState)
        Do ilStat = 1, nState
          iloc = ilStat+nState*(ilStat-1)
          WRK2(ilStat,1) = Work(ipSLag+iloc-1)
          If (.not.isNAC) Then
            Do jlStat = 1, ilStat-1
             WRK1(ilStat,jlStat)=WRK1(ilStat,jlStat)+WRK1(jlStat,ilStat)
              WRK1(jlStat,ilSTat) = 0.0d+00
            End Do
          End If
        End Do
        Call DCopy_(nSLag,[0.0D+00],0,Work(ipSLag),1)
        Call DaXpY_(nState*nState,1.0D+00,WRK1,1,Work(ipSLag),1)
      End If
C
      IF (IFXMS.or.(IFMSCOUP.and.isNAC.and.isCSF)) Then
        If (.not.IFXMS) Then
          !! For MS-CASPT2, only the second term in eq.(68)
          Call DCopy_(nState**2,[0.0D+00],0,U0,1)
          Call DCopy_(nState,[1.0D+00],0,U0,nState+1)
        End If
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        CALL XMS_Grad(Work(ipCLagFull),H0,U0,UEFF,WRK2)
        CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT =CPTF10-CPTF0
        WALLT=TIOTF10-TIOTF0
        write(6,*) "XMS_Grad: CPU/WALL TIME=", cput,wallt
      End If
C
      !! Now, compute the state Lagrangian and do some projections
      Call CLagFinal(Work(ipCLagFull),Work(ipSLag))
C
      !! Add MS-CASPT2 contributions
      If (IFMSCOUP) Then
        Do ilStat = 1, nState
          Do jlStat = 1, ilStat
            iloc = ilStat + nState*(jlStat-1)
            If (isNAC) Then
              If (.not.IFXMS.and.ilstat.ne.jlstat) Cycle
              Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
     *             + UEFF(jlStat,iRoot1)*UEFF(ilStat,iRoot2)
              Scal = Scal*0.5D+00
              Work(ipSLag+iloc-1) = Work(ipSLag+iloc-1) + Scal
C
              If (ilStat.ne.jlStat) Then
                jloc = jlStat + nState*(ilStat-1)
                Work(ipSLag+jloc-1) = Work(ipSLag+jloc-1) + Scal
              End If
            Else
              IF (IFXMS) Then
                Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                If (ilStat.ne.jlStat) Scal = Scal*2.0d+00
                Work(ipSLag+iloc-1) = Work(ipSLag+iloc-1) + Scal
              Else
                If (ilStat.eq.jlStat) Then
                  Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                  Work(ipSLag+iloc-1) = Work(ipSLag+iloc-1) + Scal
                End If
              End If
            End If
          End Do
        End Do
      End If
C
      !! Add XDW-CASPT2 contributions: derivative of <alpha|H|beta> etc
      If (IFDW) Then
        !! OMGDER is computed in XMS basis, so back-transform to CASSCF
        call sqprt(work(ipomgder),nstate)
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,U0,nState,Work(ipOMGDER),nState,
     *              0.0D+00,UEFF,nState)
        Call DGEMM_('N','T',nState,nState,nState,
     *              1.0D+00,UEFF,nState,U0,nState,
     *              0.0D+00,Work(ipOMGDER),nState)
      End If
C
      !! Subtract the original rhs_sa.f or rhs_nac.f contribution
      !! For MS-type CASPT2, CASSCF part has to be determined by UEFF
      If (IFMSCOUP.and.iRoot1.eq.iRoot2) Then
        iloc = MAX(iRoot1,iRoot2)*(MAX(iRoot1,iRoot2)-1)/2
     *       + MIN(iRoot1,iRoot2) !! iRlxRoot*(iRlxRoot+1)/2
        iloc = MAX(iRoot1,iRoot2)+nState*(MIN(iRoot1,iRoot2)-1)
        Work(ipSLag+iloc-1) = Work(ipSLag+iloc-1) - 1.0D+00
      End If
C
      If (isNAC) Then
        If (isCSF) Then
          Call CnstAntiC(Work(ipDPT2Canti),UEFF,U0)
        Else
          !! Clear just in case
          Call DCopy_(nBasSq,[0.0D+00],0,Work(ipDPT2Canti),1)
        End If
      End If
C
      !! Back-transform the CI Lagrangian
      !! It is in the XMS basis, so it has to be transformed to
      !! CASSCF basis to be used in Z-vector
      !! No need to do this for SLag.
      If (IFXMS) Then
        Call GetMem('CI1','ALLO','REAL',LCI1,nConf*nState)
        Call DGEMM_('N','T',nConf,nState,nState,
     &              1.0D+00,Work(ipCLagFull),nConf,U0,nState,
     &              0.0D+00,Work(LCI1),nConf)
        Call DCopy_(nConf*nState,Work(LCI1),1,Work(ipCLagFull),1)
        Call GetMem('CI1','FREE','REAL',LCI1,nConf*nState)
      End If
C
      DEBUG = .FALSE.
      Call Molcas_Open(LuPT2,'PT2_Lag')
C     Write (LuPT2,*) BSHIFT
      !! configuration Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "CLag"
      Do i = 1, nCLag
C       if (abs(work(ipclagfull+i-1)).le.1.0d-10) work(ipclagfull+i-1)=0.0d+00
        Write (LuPT2,*) Work(ipCLagFull+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipCLagFull+i-1)
      End Do
      !! orbital Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "OLag"
      Do i = 1, nOLag
C       if (abs(work(ipolagfull+i-1)).le.1.0d-10) work(ipolagfull+i-1)=0.0d+00
        Write (LuPT2,*) Work(ipOLagFull+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipOLagFull+i-1)
      End Do
      !! state Lagrangian (read in RHS_PT2)
      If (DEBUG) write(6,*) "SLag"
      Do i = 1, nSLag
        Write (LuPT2,*) Work(ipSLag+i-1)
        If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipSLag+i-1)
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
      If (isNAC) Then
        Do i = 1, nBasSq
C         write(6,*) i,work(ipdpt2c+i-1)
          Write (LuPT2,*) Work(ipDPT2Canti+i-1)
          If (DEBUG) write(6,'(I6,F20.10)') i,Work(ipDPT2Canti+i-1)
        End Do
      End If
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
      CALL GETMEM('CLAG   ','FREE','REAL',ipCLag   ,nCLag*2)
      CALL GETMEM('OLAG   ','FREE','REAL',ipOLag   ,nOLag*2)
      CALL GETMEM('SLAG   ','FREE','REAL',ipSLag   ,nSLag)
      CALL GETMEM('WLAG   ','FREE','REAL',ipWLag   ,nWLag)
C
      CALL GETMEM('FIFA   ','FREE','REAL',ipFIFA   ,nBasSq)
      CALL GETMEM('FIMO   ','FREE','REAL',ipFIMO   ,nBasSq)
C
      If (IFXMS) Then
        CALL GETMEM('FIFASA ','FREE','REAL',ipFIFASA ,nBasSq)
      End If
C
      If (IFDW) Then
        CALL GETMEM('OMGDER ','FREE','REAL',ipOMGDER ,nState**2)
      End If
C
      If (isNAC) Then
        Call GETMEM('DPT2Canti','FREE','REAL',ipDPT2Canti,nBasSq)
      End If
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
C
C-----------------------------------------------------------------------
C
      Subroutine GradStart
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "caspt2_grad.fh"
C
      !! CASPT2 is invariant with respect to rotations in active?
      INVAR=.TRUE.
      If (BSHIFT.NE.0.0d+00) INVAR=.FALSE.
      IF (.not.IFINVAR)      INVAR=.FALSE.
C
      If (.not.INVAR .and. IPRGLB.GE.USUAL) Then
        Write (6,*)
        Write (6,'(3X,"This is a non-invariant CASPT2 calculation")')
        If (BSHIFT.NE.0.0D+00)
     *    Write (6,'(3X,"- IPEA shift is employed")')
        Write (6,'(3X,"A linear equation will be solved to obtain ",
     *                "off-diagonal active density")')
        Write (6,*)
      End If
C
      IF (.not.IFSADREF.and.nState.ne.1.and..not.IFXMS) Then
        write(6,*)
     *    "Please add SADREF keyword in CASPT2 section",
     *    "This keyword is recommended with state-averaged reference"
C       call abend
      End If
      IF (.not.IFDORTHO.and.BSHIFT.ne.0.0D+00) Then
        write(6,*)
     *    "It seems that DORT keyword is not used, ",
     *    "even though this calculation uses the IPEA shift"
        write(6,*)
     *    "Sometimes, analytic gradients do not agree ",
     *    "with numerical gradients"
        write(6,*)
     *    "(which are correct?)"
      End If
C
      End Subroutine GradStart
C
C-----------------------------------------------------------------------
C
      Subroutine GradPrep(UEFF,VECROT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "caspt2_grad.fh"
C
C#include "nadc.fh"
C#include "nac.fh"
C
      Dimension UEFF(nState,nState),VECROT(nStLag)
C
C     H_{IJ} = <I|H^2|J>
C     we diagonalize tH_{IJ} = tilde-H_{IJ} = (H_{IJ}+H_{JI})/2
C     U^T*H*U = U_{IK}*tH_{IJ}*U_{JL}
C             = U_{IK}*(H_{IJ}+H_{JI})*U_{JL}/2 = delta_{KL}
C     Derivative of the diagonal:
C       d(U_{IK}*(H_{IJ}+H_{JI})*U_{JK})/dx/2
C       = U_{IK}*d(H_{IJ}+H_{JI})/dx*U_{JK}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JK}/2 + U_{JK}*dH_{IJ}/dx*U_{IK}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JK}
C     Derivative of H for off-diagonal:
C       d(U_{IK}*(H_{IJ}+H_{JI})*U_{JL})/dx/2
C       = U_{IK}*d(H_{IJ}+H_{JI})/dx*U_{JL}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JL}/2 + U_{JK}*dH_{IJ}/dx*U_{IL}/2
C       = (U_{IK}*U_{JL} + U_{JK}*U_{IL}) * dH_{IJ}/dx/2
C       also
C       d(U_{IL}*(H_{IJ}+H_{JI})*U_{JK})/dx/2
C       = U_{IL}*d(H_{IJ}+H_{JI})/dx*U_{JK}/2
C       = U_{IL}*dH_{IJ}/dx*U_{JK}/2 + U_{JL}*dH_{IJ}/dx*U_{IK}/2
C
C Hij = U1i*tH11*U1j + U1i*tH12*U2j + U2i*tH21*U1j + U2i*tH22*U2j
C     = (U1i*(H11+H11)*U1j + U1i*(H12+H21)*U2j + U2i*(H21+H12)*U1j + U2i*(H22+H22)*U2j)/2
C     = (2*U1i*H1j*U12 + (U1i*U2j+U2i*U1j)H12 + (U1i*U2j+U2i*U1j)*H21 + 2*U2i*H22*U2j)/2
C     If i  = j, UIi*dHij/dx*UJj
C     If i \= j, (UIi*UJj+UJi*UIj)*dHij/dx*0.5
C
      !! Construct the rotation vector
      If (IFMSCOUP) Then
        Do iState = 1, nState
          TMP = UEFF(iState,iRoot1)*UEFF(jState,iRoot2)
     *        + UEFF(iState,iRoot2)*UEFF(jState,iRoot1)
          VECROT(iState) = TMP*0.5d+00
        End Do
        jStLag    = jState
      Else
        VECROT(1) = 1.0D+00
        jStLag    = 1
      End If
C
      End Subroutine GradPrep
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFinal(OLag,Trf)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Dimension OLag(*),Trf(*)
C
      CALL GETMEM('WRK    ','ALLO','REAL',ipWRK    ,nBasSq)
      CALL GETMEM('WLAGL  ','ALLO','REAL',ipWLagL  ,nWLag)
C
      Call DCopy_(nWLag,[0.0d+00],0,Work(ipWLagL),1)
      Call DaXpY_(nWLag,0.5D+00,OLag,1,Work(ipWLagL),1)
C     write(6,*) "Wlag square"
C     call sqprt(work(ipwlag),nbast)
C
      !! W(MO) -> W(AO) using the quasi-canonical orbitals
      !! No need to back transform to natural orbital basis
      Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *            1.0D+00,Work(LCMOPT2),nBasT,Work(ipWLagL),nBasT,
     *            0.0D+00,Work(ipWRK),nBasT)
      Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *            1.0D+00,Work(ipWRK),nBasT,Work(LCMOPT2),nBasT,
     *            0.0D+00,Work(ipWLagL),nBasT)
C
      !! square -> triangle for WLag(AO)
      Call DCopy_(nBasT*nBasT,Work(ipWLagL),1,Work(ipWRK),1)
      iBasTr = 1
      iBasSq = 1
      Do iSym = 1, nSym
        nBasI = nBas(iSym)
        liBasTr = iBasTr
        liBasSq = iBasSq
        Do iBasI = 1, nBasI
          Do jBasI = 1, iBasI
            liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
            If (iBasI.eq.jBasI) Then
              Work(ipWLagL  +liBasTr-1) = Work(ipWRK  +liBasSq-1)
            Else
            liBasSq2 = iBasSq + jBasI-1 + nBasI*(iBasI-1)
              Work(ipWLagL  +liBasTr-1)
     *          = Work(ipWRK  +liBasSq-1)
     *          + Work(ipWRK  +liBasSq2-1)
            End If
            liBasTr = liBasTr + 1
          End Do
        End Do
        iBasTr = iBasTr + nBasI*(nBasI+1)/2
        iBasSq = iBasSq + nBasI*nBasI
      End Do
      Call DaXpY_(nBasSq,1.0D+00,Work(ipWLagL),1,Work(ipWLag),1)
      CALL GETMEM('WLAGL  ','FREE','REAL',ipWLagL  ,nWLag)
C
C
C
      !! Transform quasi-canonical -> natural MO basis
      !! orbital Lagrangian
      Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *            1.0D+00,Trf,nBasT,OLag,nBasT,
     *            0.0D+00,Work(ipWRK),nBasT)
      Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *            1.0D+00,Work(ipWRK),nBasT,Trf,nBasT,
     *            0.0D+00,OLag,nBasT)
      !! sufficient only for active
      nBasI = nBas(1)
      Call DCopy_(nBasI**2,OLag,1,Work(ipWRK),1)
      Call DGeSub(Work(ipWRK),nBas(1),'N',
     &            Work(ipWRK),nBas(1),'T',
     &            OLag,nBas(1),
     &            nBas(1),nBas(1))
      Call DaXpY_(nOLag,1.0D+00,OLag,1,Work(ipOLagFull),1)
C
      CALL GETMEM('WRK    ','FREE','REAL',ipWRK    ,nBasSq)
C
      End Subroutine OLagFinal
C
C-----------------------------------------------------------------------
C
      Subroutine CnstFIFAFIMO(MODE)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
C     nOrbI = nOrb(1)
C     Call SQUARE(Work(LFIFA),Work(ipFIFASA),1,nOrbI,nOrbI)
C     return
C
      If (IfChol) Then
        !! For DF or CD, we already have FIFA and FIMO in AO,
        !! so just do AO -> MO transformation
        CALL GETMEM('WRK1   ','ALLO','REAL',ipWRK1   ,nBasSq)
        CALL GETMEM('WRK2   ','ALLO','REAL',ipWRK2   ,nBasSq)
C
C           write (*,*) "ipfifa,ipfimo"
C           do i = 1, 10
C           write (*,*) i,work(ipfifa+i-1),work(ipfimo+i-1)
C           end do
        !! Read H_{\mu \nu}
C       IRC=-1
C       IOPT=6
C       ICOMP=1
C       ISYLBL=1
C       CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,Work(ipWRK1),ISYLBL)
C       Call DaXpY_(nBasTr,1.0D+00,Work(ipWRK1),1,Work(ipFIMO),1)
C       Call DaXpY_(nBasTr,1.0D+00,Work(ipWRK1),1,Work(ipFIFA),1)
C
C       Call DaXpY_(nBasTr,1.0D+00,Work(ipFIMO),1,Work(ipFIFA),1)
        iSQ = 0
        iTR = 0
        Do iSym = 1, nSym
          nOrbI = nOrb(iSym)
          nBasI = nBas(iSym)
          !! FIFA
          If (nFroT.eq.0) Then
            If (MODE.eq.0.and.IFDW) Then
              Call SQUARE(Work(LFIFA+iTr),Work(ipFIFASA+iSQ),
     *                    1,nBasI,nBasI)
            Else If (MODE.eq.1) Then
              Call SQUARE(Work(LFIFA+iTr),Work(ipFIFA+iSQ),
     *                    1,nBasI,nBasI)
C             write (*,*) "fifa in MO"
C             call sqprt(work(ipfifa+isq),nbasi)
            End If
          Else
            Call SQUARE(Work(ipFIFA+iTr),Work(ipWRK1),1,nBasI,nBasI)
C             write (*,*) "fifasa in AO"
C             call sqprt(work(ipwrk1),nbasi)
            If (MODE.eq.0.and.IFDW) Then
              !! with the state-average
              !! FIFASA will be natural basis
              Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFIFASA+iSQ),
     *                     Work(ipWRK1),Work(ipWRK2))
C             write (*,*) "fifasa in MO"
C             call sqprt(work(ipfifasa+isq),nbasi)
            Else If (MODE.eq.1) Then
              !! with the state-specific or dynamically weighted
              !! FIFA will be quasi-canonical basis
              Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFIFA+iSQ),
     *                     Work(ipWRK1),Work(ipWRK2))
C             write (*,*) "fifa in MO"
C             call sqprt(work(ipfifa+isq),nbasi)
            End If
          End If
C
          !! FIMO
C         If (MODE.eq.0) Then
            If (nFroT.eq.0) Then
              Call SQUARE(Work(LFIMO+iTr),Work(ipFIMO+iSQ),
     *                    1,nBasI,nBasI)
            Else
              !! FIMO will be natural basis
              Call SQUARE(Work(ipFIMO+iTr),Work(ipWRK1),1,nBasI,nBasI)
              Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFIMO+iSQ),
     *                     Work(ipWRK1),Work(ipOLag))
C             write (*,*) "fimo in MO"
C             call sqprt(work(ipfimo+isq),nbasi)
            End If
C         End If
          iSQ = iSQ + nBasI*nBasI
          iTR = iTR + nBasI*(nBasI+1)/2
        End Do
        CALL GETMEM('WRK1   ','FREE','REAL',ipWRK1   ,nBasSq)
        CALL GETMEM('WRK2   ','FREE','REAL',ipWRK2   ,nBasSq)
C
        If (IFXMS.and..not.IFDW)
     *    Call DCopy_(nBasSq,Work(ipFIFA),1,Work(ipFIFASA),1)
      Else
        If (nFroT.ne.0) Then
        Else
          iSQ = 0
          iTR = 0
          Do iSym = 1, nSym
            nOrbI = nOrb(iSym)
            nBasI = nBas(iSym)
            Call SQUARE(Work(LFIFA+iTr),Work(ipFIFA+iSQ),
     *                  1,nBasI,nBasI)
            Call SQUARE(Work(LFIMO+iTr),Work(ipFIMO+iSQ),
     *                  1,nBasI,nBasI)
            iSQ = iSQ + nBasI*nBasI
            iTR = iTR + nBasI*(nBasI+1)/2
          End Do
          If (IFXMS.and..not.IFDW)
     *      Call DCopy_(nBasSq,Work(ipFIFA),1,Work(ipFIFASA),1)
        End If
C
      !! XDW or RMS case: call after XDWINI
      If (MODE.eq.0) Then
      End If
C
      !! XMS case: call after GRPINI
      If (MODE.eq.1) Then
      End If
C
C     !! SS or MS case: call in dens.f
C     If (MODE.eq.2) Then
C       If (IFSADREF) Then
C       Else
C       End If
C     End If
      End If
C
      Return
C
      End Subroutine CnstFIFAFIMO
