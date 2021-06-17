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
* Copyright (C) 1998, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE CASPT2(IRETURN)
      USE SUPERINDEX
      USE INPUTDATA
      USE PT2WFN
      IMPLICIT NONE
      INTEGER IRETURN
*----------------------------------------------------------------------*
*     1998  PER-AAKE MALMQUIST                                         *
*     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
*     UNIVERSITY OF LUND, SWEDEN                                       *
*----------------------------------------------------------------------*

C     SECOND ORDER PERTURBATION CALCULATIONS WITH A CASSCF
C     REFERENCE FUNCTION.

C     ORIGINAL CASPT2 PROGRAM WRITTEN 890526 BY:
C     KERSTIN ANDERSSON (ALL CASPT2 CODE WITH FOLLOWING EXCEPTIONS)
C     PER-AKE MALMQVIST (THE GUGA PART)
C     BJORN O ROOS      (THE INTEGRAL TRANSFORMATION)

C     MODIFIED 1991-02-23 BY PER-AAKE MALMQUIST, FOR USE WITH THE
C     MOLCAS RASSCF PROGRAM.

C     ALL PROGRAM REWRITTEN (MALMQVIST 1993) FOR MOLCAS VERSION 3,
C     EXCEPT TRACTL AND TRA2 KEPT BUT SLIGHTLY MODIFIED, AND ALSO:
C     INTERFACING WITH MOLCAS-3 BY MARCUS FUELSCHER

C     REWRITTEN FOR 1) EXACT ACTIVE DENSITY MATRIX
C     2) QUANTITIES NEEDED FOR GRADIENTS
C     3) MULTI-STATE CASPT2
C     BY MALMQVIST 1998.
C
C     SINCE THEN, THE FOLLOWING MODIFICATIONS HAVE BEEN MADE:
C     IPEA HAMILTONIAN BY G. GHIGO & P. MALMQVIST
C       Chem. Phys. Lett. 396, pp 142 (2004)
C       http://dx.doi.org/10.1016/j.cplett.2004.08.032
C     LOVCASPT2/FNOCASPT2/AFREEZE BY B. ROOS & F. AQUILANTE
C       J. Chem. Phys. 131, 034113 (2009)
C       http://dx.doi.org/10.1063/1.3157463
C     CHOLESKY SUPPORT BY P. MALMQVIST AND F. AQUILANTE
C       J. Chem. Theory Comput. 4 (5), pp 694-702 (2008)
C       http://dx/doi.org/10.1021/ct700263h
C     RASPT2 BY P. MALMQVIST
C       J. Chem. Phys. 128, 204109 (2008)
C       http://dx.doi.org/10.1063/1.2920188
C     PARALLELIZATION BY S. VANCOILLIE
C       J. Comp. Chem. 34, pp 1937-1948 (2013)
C       http://dx.doi.org/10.1002/jcc.23342
C     MAJOR REFACTORING OF INITIALIZATION PHASE AND ADDITION OF
C     HDF5 SUPPORT BY S. VANCOILLIE (2016)
C
************************************************************************
#include "rasdim.fh"
#include "warnings.fh"
#include "constants.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "stdalloc.fh"
#include "caspt2_grad.fh"
      CHARACTER(60) STLNE2
#ifdef _MOLCAS_MPP_
      LOGICAL KING, Is_Real_Par
#endif
* Timers
      REAL*8  CPTF0, CPTF10, CPTF11, CPTF12, CPTF13, CPTF14,
     &       TIOTF0,TIOTF10,TIOTF11,TIOTF12,TIOTF13,TIOTF14,
     &          CPE,CPUTOT,TIOE,TIOTOT
* Indices
      INTEGER I
      INTEGER ISTATE
      INTEGER IGROUP,JSTATE_OFF
* Convergence check
      INTEGER ICONV
* Relative energies
      REAL*8  RELAU,RELEV,RELCM,RELKJ

* Effective Hamiltonian
      REAL*8, ALLOCATABLE :: Heff(:,:), Ueff(:,:)

* Zeroth-order Hamiltonian
      REAL*8, ALLOCATABLE :: H0(:,:), U0(:,:)


      !! for debug
      Logical,parameter ::  DerMO=.false.,
     *                      DerCI=.false.,
     *                      DerGen=.false.,
     *                      DerICB=.false.
      Integer iOrb,jOrb,iAO,iMO,iVib,ipWRK1,ipWRK2,ipWRK3,indCI,idCI,
     *        idT,iRoots,jRoots,ijRoots,ipFIMO,ipFIFA,nAS,nIN,icase,isym
      Real*8  Delta,EForward,EBackward,EigI,EigJ,SCAL,SLagV
      Real*8, External :: DDot_

      Call StatusLine('CASPT2:','Just starting')

      IRETURN = 0
      CALL QENTER('CASPT2')

      CALL SETTIM
      ! CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)

* Probe the environment to globally set the IPRGLB value
      Call Set_Print_Level

*======================================================================*
*
      Call StatusLine('CASPT2:','Initializing')
      CALL PT2INI
* Initialize effective Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(Heff,Nstate,Nstate)
      CALL MMA_ALLOCATE(Ueff,Nstate,Nstate)
      Heff=0.0D0
      Ueff=0.0D0
* Initialize zeroth-order Hamiltonian and eigenvectors
      CALL MMA_ALLOCATE(H0,Nstate,Nstate)
      CALL MMA_ALLOCATE(U0,Nstate,Nstate)
      H0=0.0D0
* U0 is initialized as the identity matrix, in the case of a
* standard MS-CASPT2 calculation it will not be touched anymore
      U0=0.0D0
      call dcopy_(Nstate,[1.0d0],0,U0,Nstate+1)
*
*======================================================================*
* If the EFFE keyword has been used, we already have the multi state
* coupling Hamiltonian effective matrix, just copy the energies and
* proceed to the MS coupling section.
* Otherwise, put the CASSCF energies on the diagonal, i.e. form the
* first-order corrected effective Hamiltonian:
*     Heff[1] = PHP
* and later on we will add the second-order correction
* Heff(2) = PH \Omega_1 P to Heff[1]
      IF (INPUT%JMS) THEN
        DO I=1,NSTATE
          ENERGY(I)=INPUT%HEFF(I,I)
        END DO
        HEFF(:,:)=INPUT%HEFF(:,:)
        GOTO 1000
      ELSE
        DO I=1,NSTATE
          HEFF(I,I) = REFENE(I)
        END DO
        IF (IPRGLB.GE.DEBUG) THEN
          write(6,*)' Heff[1] in the original model space basis:'
          call prettyprint(Heff,Nstate,Nstate)
        END IF
      END IF

* In case of a XDW-CASPT2 calculation we first rotate the CASSCF
* states according to the XMS prescription in xdwinit
      if (IFXMS.and.IFDW) then
        call xdwinit(Heff,H0,U0)
        if (IFEFOCK) then
          call wgtinit(H0)
        else
          call wgtinit(Heff)
        end if
      else
        call wgtinit(Heff)
      end if

* Before entering the long loop over groups and states, precompute
* the 1-RDMs for all states and mix them according to the type of
* calculation: MS, XMS, DW or XDW.
      call rdminit

* For (X)Multi-State, a long loop over root states.
* The states are ordered by group, with each group containing a number
* of group states for which GRPINI is called.
      JSTATE_OFF=0
      STATELOOP: DO IGROUP=1,NGROUP

       IF ((NLYGROUP.NE.0).AND.(IGROUP.NE.NLYGROUP)) THEN
         JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
         CYCLE
       END IF

       IF (IPRGLB.GE.USUAL) THEN
        If(.not.IFNOPT2) Then
         WRITE(STLNE2,'(A,1X,I3)') 'CASPT2 computation for group',IGROUP
         CALL CollapseOutput(1,TRIM(STLNE2))
         WRITE(6,*)
        End If
       END IF

       CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
C      write(6,*) "before grpini"
       CALL GRPINI(IGROUP,NGROUPSTATE(IGROUP),JSTATE_OFF,HEFF,H0,U0)
C      write(6,*) "after grpini"
       CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
       CPUGIN=CPTF10-CPTF0
       TIOGIN=TIOTF10-TIOTF0
CProducing XMS Rotated States
       If(IFNOPT2) then
        GOTO 9999
       END IF

       DO ISTATE=1,NGROUPSTATE(IGROUP)
         JSTATE = JSTATE_OFF + ISTATE


* Skip this state if we only need 1 state and it isn't this "one".
         IF ((NLYROOT.NE.0).AND.(JSTATE.NE.NLYROOT)) CYCLE

         CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
C      write(6,*) "before stini"
C            CALL GETMEM('WRK1','ALLO','REAL',ipWRK1,nConf*nRoots)
C            idCI = idTCEX
C            Call dDaFile(LUCIEX,2,Work(ipWRK1),nConf,idCI)
C            Work(ipWRK1) = Work(ipWRK1) + 1.0D-05
C            idCI = idTCEX
C            Call dDaFile(LUCIEX,1,Work(ipWRK1),nConf,idCI)
C            CALL GETMEM('WRK1','FREE','REAL',ipWRK1,nConf*nRoots)
         CALL STINI
         CALL TIMING(CPTF11,CPE,TIOTF11,TIOE)
         CPUSIN=CPTF11-CPTF0
         TIOSIN=TIOTF11-TIOTF0

* Solve CASPT2 equation system and compute corr energies.
         IF (IPRGLB.GE.USUAL) THEN
            WRITE(6,'(20A4)')('****',I=1,20)
            WRITE(6,*)' CASPT2 EQUATION SOLUTION'
            WRITE(6,'(20A4)')('----',I=1,20)
         END IF

         Write(STLNE2,'(A27,I3)')'Solve CASPT2 eqs for state ',
     &                               MSTATE(JSTATE)
         Call StatusLine('CASPT2:',TRIM(STLNE2))
         CALL EQCTL2(ICONV)
C
C     if (.not.ifdens) then
C     write(6,*) "quit for numerical derivatives"
C     call abend()
C     end if

        DoPT2NUM = .false.
         If (DerMO.and.iRLXroot.eq.jstate) Then
           IPRGLB = SILENT
           Delta = 1.0d-05
           write(6,*) "Numerical Orbital Lagrangian"
           write(6,*) "norbt,nbast = ", norbt,nbast
           write(6,*) "Delta = ", Delta
           CALL GETMEM('WRK1','ALLO','REAL',ipWRK1,nBasT*nBasT)
           CALL GETMEM('WRK2','ALLO','REAL',ipWRK2,nBasT*nBasT)
           CALL GETMEM('LCMO','ALLO','REAL',LCMO,nBasT*nBasT)
           CALL GETMEM('LFMO','ALLO','REAL',ipFIMO,NFIMO)
           CALL GETMEM('LFMO','ALLO','REAL',ipFIFA,NFIFA)
           Call DCopy_(nBasT*nBasT,[0.0D+00],0,Work(ipWRK1),1)
           Call DCopy_(NFIMO,Work(LFIMO),1,Work(ipFIMO),1)
           Call DCopy_(NFIFA,Work(LFIFA),1,Work(ipFIFA),1)
           Do iMO = nFro(1)+1, nFro(1)+nOrbT
             Do iAO = 1, nBasT
               Do iVib = 1, 2
                 If (iVib.eq.1) Then
                   Work(LCMOPT2+iAO-1+nBast*(iMO-1))
     *               = Work(LCMOPT2+iAO-1+nBast*(iMO-1)) + Delta
                 Else If (iVib.eq.2) Then
                   Work(LCMOPT2+iAO-1+nBasT*(iMO-1))
     *               = Work(LCMOPT2+iAO-1+nBasT*(iMO-1)) - Delta
                 End If
                 Call DCopy_(nBasT*nBasT,Work(LCMOPT2),1,Work(LCMO),1)
C
                 CALL STINI
                 Call TraCtl(0)
C                Call DCopy_(NFIMO,Work(ipFIMO),1,Work(LFIMO),1)
                 Call DCopy_(NFIFA,Work(ipFIFA),1,Work(LFIFA),1)
                 Call EqCtl2(iConv)
C
                 If (iVib.eq.1) Then
                   EForward = E2Tot
                   Work(LCMOPT2+iAO-1+nBasT*(iMO-1))
     *               = Work(LCMOPT2+iAO-1+nBasT*(iMO-1)) - Delta
                 Else If (iVib.eq.2) Then
                   EBackward = E2Tot
                   Work(LCMOPT2+iAO-1+nBasT*(iMO-1))
     *               = Work(LCMOPT2+iAO-1+nBasT*(iMO-1)) + Delta
                 End If
                 Call DCopy_(nBasT*nBasT,Work(LCMOPT2),1,Work(LCMO),1)
               End Do
               Work(ipWRK1+iAO-1+nBasT*(iMO-1))
     *           = (EForward-EBackward)/(2.0D+00*Delta)
             End Do
           End Do
C
           Call DGEMM_('T','N',nFro(1)+nOrbT,nFro(1)+nOrbT,nBasT,
     *                 1.0D+00,Work(LCMOPT2),nBasT,Work(ipWRK1),nBasT,
     *                 0.0D+00,Work(ipWRK2),nFro(1)+nOrbT)
           write(6,*) "Orbital Lagrangian"
           Call SqPrt(Work(ipWRK2),nFro(1)+nOrbT)
C
           write(6,*) "Anti-Symmetrized Orbital Lagrangian"
           Call DGeSub(Work(ipWRK2),nFro(1)+nOrbT,'N',
     *                 Work(ipWRK2),nFro(1)+nOrbT,'T',
     *                 Work(ipWRK1),nFro(1)+nOrbT,
     *                 nFro(1)+nOrbT,nFro(1)+nOrbT)
           Call SqPrt(Work(ipWRK1),nFro(1)+nOrbT)
C
           write(6,*) "Anti-Symmetrized Orbital Lagrangian/Eigenvalues"
           write(6,*) "This should be DPT2 in MO"
           Do iOrb = 1, nFro(1)+nOrbT
             If (iOrb.le.nFro(1)) Then
               call getfdiag(iorb,eigi,nfro(1)+norbt)
             Else
               EigI = EPS(iOrb-nFro(1))
             End If
             write(6,'("EPS(",i2,") = ",f20.10)') iorb,eigi
             Do jOrb = 1, iOrb-1
               If (jOrb.le.nFro(1)) Then
                 call getfdiag(jorb,eigj,nfro(1)+norbt)
               Else
                 EigJ = EPS(jOrb-nFro(1))
               End If
               Work(ipWRK1+iOrb-1+(nFro(1)+nOrbT)*(jOrb-1))
     *        = Work(ipWRK1+iOrb-1+(nFro(1)+nOrbT)*(jOrb-1))/(EigJ-EigI)
               Work(ipWRK1+jOrb-1+(nFro(1)+nOrbT)*(iOrb-1))
     *        =-Work(ipWRK1+jOrb-1+(nFro(1)+nOrbT)*(iOrb-1))/(EigJ-EigI)
             End Do
           End Do
           Call DScal_((nFro(1)+nOrbT)**2,0.5D+00,Work(ipWRK1),1)
           !! Orbital energy derivative
           Call TraCtl(0)
           Do iMO = 1, nOrbT
C            If (iMO.gt.nIsh(1).and.iMO.le.nIsh(1)+nAsh(1)) Cycle
             Do iVib = 1, 2
               If (iVib.eq.1) Then
                 EPS(iMO) = EPS(iMO) + Delta
                 If (iMO.le.nIsh(1)) Then
                   EPSI(iMO) = EPSI(iMO) + Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)) Then
                   EPSA(iMO-nIsh(1)) = EPSA(iMO-nIsh(1)) + Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)+nSsh(1)) Then
                   EPSE(iMO-nIsh(1)-nAsh(1))
     *               = EPSE(iMO-nIsh(1)-nAsh(1)) + Delta
                 End If
               Else If (iVib.eq.2) Then
                 EPS(iMO) = EPS(iMO) - Delta
                 If (iMO.le.nIsh(1)) Then
                   EPSI(iMO) = EPSI(iMO) - Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)) Then
                   EPSA(iMO-nIsh(1)) = EPSA(iMO-nIsh(1)) - Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)+nSsh(1)) Then
                   EPSE(iMO-nIsh(1)-nAsh(1))
     *               = EPSE(iMO-nIsh(1)-nAsh(1)) - Delta
                 End If
               End If
C
               CALL STINI
               Call EqCtl2(iConv)
C              if (imo.eq.12) then
C                write(6,*) "asdf"
C                write(6,*) e2tot
C              end if
C                write(6,*) "nin = ",nindep(1,4)
C
               If (iVib.eq.1) Then
                 EForward = E2Tot
                 EPS(iMO) = EPS(iMO) - Delta
                 If (iMO.le.nIsh(1)) Then
                   EPSI(iMO) = EPSI(iMO) - Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)) Then
                   EPSA(iMO-nIsh(1)) = EPSA(iMO-nIsh(1)) - Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)+nSsh(1)) Then
                   EPSE(iMO-nIsh(1)-nAsh(1))
     *               = EPSE(iMO-nIsh(1)-nAsh(1)) - Delta
                 End If
               Else If (iVib.eq.2) Then
                 EBackward = E2Tot
                 EPS(iMO) = EPS(iMO) + Delta
                 If (iMO.le.nIsh(1)) Then
                   EPSI(iMO) = EPSI(iMO) + Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)) Then
                   EPSA(iMO-nIsh(1)) = EPSA(iMO-nIsh(1)) + Delta
                 Else If (iMO.le.nIsh(1)+nAsh(1)+nSsh(1)) Then
                   EPSE(iMO-nIsh(1)-nAsh(1))
     *               = EPSE(iMO-nIsh(1)-nAsh(1)) + Delta
                 End If
               End If
             End Do
             Work(ipWRK1+nFro(1)+iMO-1+nBasT*(nFro(1)+iMO-1))
     *         = (EForward-EBackward)/(2.0D+00*Delta)
           End Do
           Call SqPrt(Work(ipWRK1),nFro(1)+nOrbT)
C
           write(6,*) "ipWRK1,ipWRK2,nBasT"
           write(6,*) ipWRK1,ipWRK2,nBasT
           CALL GETMEM('WRK1','FREE','REAL',ipWRK1,nBasT*nBasT)
           CALL GETMEM('WRK2','FREE','REAL',ipWRK2,nBasT*nBasT)
           CALL GETMEM('LCMO','FREE','REAL',LCMO,nBasT*nBasT)
           CALL GETMEM('LFMO','FREE','REAL',ipFIMO,NFIMO)
           CALL GETMEM('LFMO','FREE','REAL',ipFIFA,NFIFA)
           Call AbEnd()
         End If





         If (DerCI.and.iRLXroot.eq.jstate) Then
           IPRGLB = SILENT
           Delta = 1.0d-06
           write(6,*) "Numerical Configuration Lagrangian"
           write(6,*) "nConf, nRoots = ", nConf,nRoots
           write(6,*) "Delta = ", Delta
           CALL GETMEM('WRK1','ALLO','REAL',ipWRK1,nConf*nRoots)
           CALL GETMEM('WRK2','ALLO','REAL',ipWRK2,nConf*nRoots)
           CALL GETMEM('LCMO','ALLO','REAL',LCMO,nBasT*nBasT)
           If (ISCF.EQ.0) Then
             idCI = idTCEX
             Call dDaFile(LUCIEX,2,Work(ipWRK1),nConf*nRoots,idCI)
           Else
             Work(ipWRK1)=1.0D+00
           End If
           Call DCopy_(nBasT*nBasT,Work(LCMOPT2),1,Work(LCMO),1)
           Do indCI = 1, Min(nConf*nRoots,500)
             Do iVib = 1, 2
C              Call GRPINI
               If (iVib.eq.1) Then
                 Work(ipWRK1+indCI-1) = Work(ipWRK1+indCI-1) + Delta
               Else If (iVib.eq.2) Then
                 Work(ipWRK1+indCI-1) = Work(ipWRK1+indCI-1) - Delta
               End If
               idCI = idTCEX
               Call dDaFile(LUCIEX,1,Work(ipWRK1),nConf*nRoots,idCI)
C
               CALL STINI
C              Call TraCtl(0)
               Call EqCtl2(iConv)
C              write(6,*) e2tot
C              write(6,*) "nin = ",nindep(1,4)
C
               If (iVib.eq.1) Then
                 EForward  = E2Tot
                 Work(ipWRK1+indCI-1) = Work(ipWRK1+indCI-1) - Delta
               Else If (iVib.eq.2) Then
                 EBackward = E2Tot
                 Work(ipWRK1+indCI-1) = Work(ipWRK1+indCI-1) + Delta
               End If
C              Call dDaFile(LUCIEX,1,Work(ipWRK1),nConf*nRoots,idCI)
             End DO
             Work(ipWRK2+indCI-1) = (EForward-EBackward)/(2.0D+00*Delta)
             write(6,'(i5,f20.10)') indCI,Work(ipWRK2+indCI-1)
           End Do
C
C          write(6,*) "Configuration Lagrangian"
C          Do indCI = 1, nConf*nRoots
C            write(6,'(i5,f20.10)') indCI,Work(ipWRK2+indCI-1)
C          End Do
           If (nRoots.ne.1) Then
             write(6,*) "Numerical State Lagrangian when nRoots > 1"
             ijRoots = 0
             Do iRoots = 1, nRoots
               Do jRoots = 1, iRoots-1
                 ijRoots = ijRoots + 1
                 SLagV = DDot_(nConf,Work(ipWRK1+nConf*(iRoots-1)),1,
     *                               Work(ipWRK2+nConf*(jRoots-1)),1)
     *                 - DDot_(nConf,Work(ipWRK1+nConf*(jRoots-1)),1,
     *                               Work(ipWRK2+nConf*(iRoots-1)),1)
                 write(6,'(2i3,f20.10)') iRoots,jRoots,SLagV
               End Do
             End Do
           End If
C          write(6,*) "After Projection"
C          Do iRoots = 1, nRoots
C            write(6,*) "iRoots = ", iRoots
C            Do jRoots = 1, nRoots
C              Scal = DDot_(nConf,Work(ipWRK1+nConf*(jRoots-1)),1,
C    *                            Work(ipWRK2+nConf*(iRoots-1)),1)
C              Call DaXpY_(nConf,-Scal,Work(ipWRK1+nConf*(jRoots-1)),1,
C    *                                 Work(ipWRK2+nConf*(iRoots-1)),1)
C            End Do
C            Do indCI = 1, nConf
C              write(6,'(i5,f20.10)') indCI,Work(ipWRK2+indCI-1)
C            End Do
C          End Do
C
           write(6,*)
           write(6,*) "Numerical EASUM"
           If (ISCF.EQ.0) Then
             idCI = idTCEX
             Call dDaFile(LUCIEX,1,Work(ipWRK1),nConf,idCI)
           Else
             Work(ipWRK1)=1.0D+00
           End If
           CALL STINI
           !! Forward
           EASUM = EASUM + Delta
           Call EqCtl2(iConv)
           EASUM = EASUM - Delta
           EForward = E2Tot
           !! Backward
           CALL STINI
           EASUM = EASUM - Delta
           Call EqCtl2(iConv)
           EASUM = EASUM + Delta
           EBackward = E2Tot
           write(6,'("DEASUM = ",f20.10)')
     *       (EForward-EBackward)/(2.0D+00*delta)
C
           write(6,*)
           write(6,*) "Numerical DEPSA"
           Do iMO = 1, nAshT
             Do iVib = 1, 2
               If (iVib.eq.1) Then
                 EPSA(iMO) = EPSA(iMO) + Delta
               Else If (iVib.eq.2) Then
                 EPSA(iMO) = EPSA(iMO) - Delta
               End If
C
               CALL STINI
               Call EqCtl2(iConv)
C
               If (iVib.eq.1) Then
                 EForward  = E2Tot
                 EPSA(iMO) = EPSA(iMO) - Delta
               Else If (iVib.eq.2) Then
                 EBackward = E2Tot
                 EPSA(iMO) = EPSA(iMO) + Delta
               End If
             End Do
             write(6,'(i3,f20.10)') iMO,
     *         (EForward-EBackward)/(2.0D+00*Delta)
           End Do
C
           CALL GETMEM('WRK1','FREE','REAL',ipWRK1,nConf*nRoots)
           CALL GETMEM('WRK2','FREE','REAL',ipWRK2,nConf*nRoots)
           CALL GETMEM('LCMO','FREE','REAL',LCMO,nBasT*nBasT)
           Call AbEnd()
         End If



         If (DerGen) Then
           IPRGLB = SILENT
           write(6,*) "Numerical ... Something"
           write(6,*) "TITLE: diagonal B derivative"
           Delta = 1.0d-05
           PT2Delta = Delta
           DoPT2Num = .true.

           CALL GETMEM('LFMO','ALLO','REAL',ipFIMO,NFIMO)
           Call DCopy_(NFIMO,Work(LFIMO),1,Work(ipFIMO),1)

           DoTriPT2 = .true.
           Do iDiffPT2 = 1, nbast*(nbast+1)/2
             Do iVibPT2 = 1, 2
C     !! template
C     If (DoPT2Num) Then
C       If (iVibPT2.eq.1) Then
C         Work(ipTEST+iDiffPT2-1) = Work(ipTEST+iDiffPT2-1) + PT2Delta
C       Else
C         Work(ipTEST+iDiffPT2-1) = Work(ipTEST+iDiffPT2-1) - PT2Delta
C       End If
C     End If
C
C              CALL STINI
C              Call TraCtl(0)
      If (DoPT2Num) Then
        Call DCopy_(NFIMO,Work(ipFIMO),1,Work(LFIMO),1)
        If (iVibPT2.eq.1) Then
          Work(LFIMO+iDiffPT2-1) = Work(LFIMO+iDiffPT2-1) + PT2Delta
        Else
          Work(LFIMO+iDiffPT2-1) = Work(LFIMO+iDiffPT2-1) - PT2Delta
        End If
      End If
               Call EqCtl2(iConv)
C
               If (iVibPT2.eq.1) EForward  = E2Tot
               If (iVibPT2.eq.2) EBackward = E2Tot
             End DO
             SCAL = (EForward-EBackward)/(2.0D+00*Delta)
             !! i*(i+1)/2 = ind
             !! i^2 + i - 2*ind = 0
             !! (-1 \pm sqrt(4*ind*ind-1))/2
             SLagV = (sqrt(1.0d+00+dble(8*iDiffPT2))-1.0d+00)*0.5d+00
C     write(6,'(i3,3f20.10)') idiffpt2,SLagV, floor(slagv+1.0d-03),
C    *    abs(slagv-floor(slagv+1.0d-03))
             If (DoTriPT2.and.
     *         abs(SLagV-floor(SLagV+1.0d-05)).gt.1.0d-08) Then
               SCAL = SCAL*0.5d+00
               write(6,'(i3,f20.10," *")') iDiffPT2,SCAL
             Else
               write(6,'(i3,f20.10)') iDiffPT2,SCAL
             End If
           End Do
           CALL GETMEM('LFMO','FREE','REAL',ipFIMO,NFIMO)
           Call AbEnd()
         End If



         If (DerICB) Then
           IPRGLB = SILENT
           write(6,*) "Numerical ICB derivative"
           Delta = 1.0d-05
           PT2Delta = Delta
           DoPT2Num = .true.

           iCase = 1
           iSym  = 1
           nIN = nINDEP(iSym,iCase)
           nAS = nASUP(iSym,iCase)

           !! derivative
           CALL GETMEM('WRK1','ALLO','REAL',ipWRK1,nAS*nAS)
           !! trans
           CALL GETMEM('WRK2','ALLO','REAL',ipWRK2,nAS*nAS)
           !! wrk
           CALL GETMEM('WRK3','ALLO','REAL',ipWRK3,nAS*nAS)

           Call DCopy_(nAS*nAS,[0.0D+00],0,Work(ipWRK1),1)
           idT  = idTMAT(iSym,iCase)
           CALL DDAFILE(LUSBT,2,Work(ipWRK2),nAS*nIN,idT)
           iDiffPT2 = 0
           Do iMO = 1, nIN
             Do iAO = 1, nAS
               iDiffPT2 = iDiffPT2 + 1
               Do iVibPT2 = 1, 2
C     !! template
C     If (DoPT2Num) Then
C       If (iVibPT2.eq.1) Then
C         Work(ipTEST+iDiffPT2-1) = Work(ipTEST+iDiffPT2-1) + PT2Delta
C       Else
C         Work(ipTEST+iDiffPT2-1) = Work(ipTEST+iDiffPT2-1) - PT2Delta
C       End If
C     End If
C
                 Call EqCtl2(iConv)
C
                 If (iVibPT2.eq.1) EForward  = E2Tot
                 If (iVibPT2.eq.2) EBackward = E2Tot
               End Do
               Work(ipWRK1+iAO-1+nAS*(iMO-1))
     *           = (EForward-EBackward)/(2.0D+00*Delta)
             End Do
           End Do
C
           !! Restore the original eigenvalues and vectors
           Call MKBMAT
           Call SBDIAG
C
           Call DGEMM_('T','N',nIN,nIN,nAS,
     *                 1.0D+00,Work(ipWRK2),nAS,Work(ipWRK1),nAS,
     *                 0.0D+00,Work(ipWRK3),nIN)
           write(6,*) "Orbital Lagrangian in ICB"
           Call SqPrt(Work(ipWRK3),nIN)
C
           write(6,*) "Anti-Symmetrized Orbital Lagrangian"
           Call DGeSub(Work(ipWRK3),nIN,'N',
     *                 Work(ipWRK3),nIN,'T',
     *                 Work(ipWRK1),nIN,nIN,nIN)
           Call SqPrt(Work(ipWRK1),nIN)
C
           write(6,*) "Anti-Symmetrized Orbital Lagrangian/Eigenvalues"
           write(6,*) "This should be density in ICB"
           idT  = idBMAT(iSym,iCase)
           CALL DDAFILE(LUSBT,2,Work(ipWRK2),nIN,idT)
           Do iOrb = 1, nIN
             EigI = Work(ipWRK2+iOrb-1)
             write(6,'("EPS(",i2,") = ",f20.10)') iorb,eigi
             Do jOrb = 1, iOrb-1
               EigJ = Work(ipWRK2+jOrb-1)
               Work(ipWRK1+iOrb-1+nIN*(jOrb-1))
     *        = Work(ipWRK1+iOrb-1+nIN*(jOrb-1))/(EigJ-EigI)
               Work(ipWRK1+jOrb-1+nIN*(iOrb-1))
     *        =-Work(ipWRK1+jOrb-1+nIN*(iOrb-1))/(EigJ-EigI)
             End Do
           End Do
           Call DScal_(nIN**2,0.5D+00,Work(ipWRK1),1)
C
           Do iMO = 1, nIN
             Do iVib = 1, 2
               If (iVib.eq.1) Then
                 Work(ipWRK2+iMO-1) = Work(ipWRK2+iMO-1) + Delta
               Else If (iVib.eq.2) Then
                 Work(ipWRK2+iMO-1) = Work(ipWRK2+iMO-1) - Delta
               End If
C
               idT  = idBMAT(iSym,iCase)
               CALL DDAFILE(LUSBT,1,Work(ipWRK2),nIN,idT)
               CALL NADIAG
               CALL RHS_INIT
               CALL MKRHS(IVECW)
               CALL PTRTOSR(1,IVECW,IRHS)
               CALL PCG(ICONV)
C
               If (iVib.eq.1) Then
                 EForward = E2Tot
                 Work(ipWRK2+iMO-1) = Work(ipWRK2+iMO-1) - Delta
               Else If (iVib.eq.2) Then
                 EBackward = E2Tot
                 Work(ipWRK2+iMO-1) = Work(ipWRK2+iMO-1) + Delta
               End If
             End Do
             Work(ipWRK1+iMO-1+nIN*(iMO-1))
     *         = (EForward-EBackward)/(2.0D+00*Delta)
           End Do
C
           Call SqPrt(Work(ipWRK1),nIN)

           !! derivative
           CALL GETMEM('WRK1','FREE','REAL',ipWRK1,nAS*nAS)
           !! trans
           CALL GETMEM('WRK2','FREE','REAL',ipWRK2,nAS*nAS)
           !! wrk
           CALL GETMEM('WRK3','FREE','REAL',ipWRK3,nAS*nAS)
           Call AbEnd()
         End If

* Save the final caspt2 energy in the global array ENERGY():
         ENERGY(JSTATE)=E2TOT

         CALL TIMING(CPTF12,CPE,TIOTF12,TIOE)
         CPUPT2=CPTF12-CPTF11
         TIOPT2=TIOTF12-TIOTF11

         IF (ICONV .NE. 0) THEN
C     No convergence. Skip the rest of the calculation.
            IRETURN = _RC_NOT_CONVERGED_
            EXIT STATELOOP
         END IF

C     Orbitals, properties:

         ! if the dens keyword is used, need accurate density and
         ! for that the serial LUSOLV file is needed, in that case copy
         ! the distributed LURHS() to LUSOLV here.
         IF(IFDENS.and.iRLXroot.eq.MSTATE(JSTATE)) THEN
           CALL PCOLLVEC(IRHS,0)
           CALL PCOLLVEC(IVECX,0)
           CALL PCOLLVEC(IVECR,0)
           CALL PCOLLVEC(IVECC,1)
           CALL PCOLLVEC(IVECC2,1)
           CALL PCOLLVEC(IVECW,1)

           CALL GrdIni
         END IF

         IF (IFPROP.and.IRLXroot.eq.MSTATE(JSTATE)) THEN
           IF (IPRGLB.GE.USUAL) THEN
             WRITE(6,*)
             WRITE(6,'(20A4)')('****',I=1,20)
             WRITE(6,*)' CASPT2 PROPERTY SECTION'
           END IF
           CALL PRPCTL
         ELSE
           IF (IPRGLB.GE.USUAL) THEN
             WRITE(6,*)
             WRITE(6,*)'  (Skipping property calculation,'
             WRITE(6,*)'   use PROP keyword to activate)'
           END IF
         END IF

         CALL TIMING(CPTF13,CPE,TIOTF13,TIOE)
         CPUPRP=CPTF13-CPTF12
         TIOPRP=TIOTF13-TIOTF12

C     Gradients.
C     Note: Quantities computed in gradients section can also
C     be used efficiently for computing Multi-State HEFF.
C     NOTE: atm the MS-CASPT2 couplings computed here are wrong!
         IF(IFDENS.and.IRLXroot.eq.MSTATE(JSTATE)) THEN
           IF (IPRGLB.GE.VERBOSE) THEN
              WRITE(6,*)
              WRITE(6,'(20A4)')('****',I=1,20)
              IF(NSTATE.GT.1) THEN
              WRITE(6,*)' CASPT2 GRADIENT/MULTI-STATE COUPLINGS SECTION'
              ELSE
                 WRITE(6,*)' CASPT2 GRADIENT SECTION'
              END IF
           END IF
           Call StatusLine('CASPT2:','Multi-State couplings')
C-SVC: for now, this part is only performed on the master node
#ifdef _MOLCAS_MPP_
           IF (Is_Real_Par()) THEN
             Call Set_Do_Parallel(.False.)
             IF (KING()) CALL GRDCTL(HEFF)
             Call Set_Do_Parallel(.True.)
             CALL GASync
           ELSE
             CALL GRDCTL(HEFF)
           END IF
#else
           CALL GRDCTL(HEFF)
#endif
           Call GrdCls
         END IF

         IF((.NOT.IFDENS) .AND. IFMSCOUP) THEN
C     If this was NOT a gradient, calculation, then the multi-state
C     couplings are more efficiently computed via three-body
C     transition density matrices.
           IF (IPRGLB.GE.VERBOSE) THEN
              WRITE(6,*)
              WRITE(6,'(20A4)')('****',I=1,20)
              WRITE(6,*)' CASPT2 MULTI-STATE COUPLINGS SECTION'
           END IF
           Call StatusLine('CASPT2:','Multi-State couplings')
           CALL MCCTL(HEFF)
         END IF


         CALL TIMING(CPTF14,CPE,TIOTF14,TIOE)
         CPUGRD=CPTF14-CPTF13
         TIOGRD=TIOTF14-TIOTF13

         CPUTOT=CPTF14-CPTF0
         TIOTOT=TIOTF14-TIOTF0

         IF (ISTATE.EQ.1) THEN
           CPUTOT=CPUTOT+CPUGIN
           TIOTOT=TIOTOT+TIOGIN
         ELSE
           CPUGIN=0.0D0
           TIOGIN=0.0D0
           CPUFMB=0.0D0
           TIOFMB=0.0D0
           CPUINT=0.0D0
           TIOINT=0.0D0
         END IF

C       IF (IPRGLB.GE.VERBOSE) THEN
        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,'(A,I6)')    '  CASPT2 TIMING INFO FOR STATE ',
     &                         MSTATE(JSTATE)
          WRITE(6,*)
          WRITE(6,'(A)')       '                        '//
     &                         ' cpu time  (s) '//
     &                         ' wall time (s) '
          WRITE(6,'(A)')       '                        '//
     &                         ' ------------- '//
     &                         ' ------------- '
          WRITE(6,*)
          WRITE(6,'(A,2F14.2)')'  Group initialization  ',CPUGIN,TIOGIN
          WRITE(6,'(A,2F14.2)')'  - Fock matrix build   ',CPUFMB,TIOFMB
          WRITE(6,'(A,2F14.2)')'  - integral transforms ',CPUINT,TIOINT
          WRITE(6,'(A,2F14.2)')'  State initialization  ',CPUSIN,TIOSIN
          WRITE(6,'(A,2F14.2)')'  - density matrices    ',CPUFG3,TIOFG3
          WRITE(6,'(A,2F14.2)')'  CASPT2 equations      ',CPUPT2,TIOPT2
          WRITE(6,'(A,2F14.2)')'  - H0 S/B matrices     ',CPUSBM,TIOSBM
          WRITE(6,'(A,2F14.2)')'  - H0 S/B diag         ',CPUEIG,TIOEIG
          WRITE(6,'(A,2F14.2)')'  - H0 NA diag          ',CPUNAD,TIONAD
          WRITE(6,'(A,2F14.2)')'  - RHS construction    ',CPURHS,TIORHS
          WRITE(6,'(A,2F14.2)')'  - PCG solver          ',CPUPCG,TIOPCG
          WRITE(6,'(A,2F14.2)')'    - scaling           ',CPUSCA,TIOSCA
          WRITE(6,'(A,2F14.2)')'    - lin. comb.        ',CPULCS,TIOLCS
          WRITE(6,'(A,2F14.2)')'    - inner products    ',CPUOVL,TIOOVL
          WRITE(6,'(A,2F14.2)')'    - basis transforms  ',CPUVEC,TIOVEC
          WRITE(6,'(A,2F14.2)')'    - sigma routines    ',CPUSGM,TIOSGM
          WRITE(6,'(A,2F14.2)')'  - array collection    ',CPUSER,TIOSER
          WRITE(6,'(A,2F14.2)')'  Properties            ',CPUPRP,TIOPRP
          WRITE(6,'(A,2F14.2)')'  Gradient/MS coupling  ',CPUGRD,TIOGRD
          WRITE(6,'(A,2F14.2)')' Total time             ',CPUTOT,TIOTOT
          WRITE(6,*)
        END IF

C End of long loop over states in the group
       END DO
       IF (IPRGLB.GE.USUAL) THEN
        If(.not.IFNOPT2) Then
         CALL CollapseOutput(0,'CASPT2 computation for group ')
         WRITE(6,*)
        End If
       END IF
C End of long loop over groups
        JSTATE_OFF = JSTATE_OFF + NGROUPSTATE(IGROUP)
9999    write (6,*)
      END DO STATELOOP

1000  CONTINUE

      IF (IRETURN.NE.0) GOTO 9000
       If(IFNOPT2) then  !XMS Skip multistate calculation.
        write(6,*)'PT2 calculation skipped with XROH keyword'
        write(6,*)
        CALL MMA_DEALLOCATE(UEFF)
        CALL MMA_DEALLOCATE(U0)
       Else
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,*)' Total CASPT2 energies:'
       DO I=1,NSTATE
        IF ((NLYROOT.NE.0).AND.(I.NE.NLYROOT)) CYCLE
        CALL PrintResult(6,'(6x,A,I3,5X,A,F16.8)',
     &    'CASPT2 Root',MSTATE(I),'Total energy:',ENERGY(I),1)
       END DO
       WRITE(6,*)
       IF (IFXMS) THEN
        WRITE(6,*)' Note that these CASPT2 energies are obtained using'
        WRITE(6,*)' the XMS Fock operator and thus do not correspond'
        WRITE(6,*)' to the true single-state CASPT2 ones.'
       END IF
       WRITE(6,*)
      END IF
      IF(IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
       WRITE(6,*)' Relative CASPT2 energies:'
       WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &   'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
       ISTATE=1
       DO I=2,NSTATE
         IF (ENERGY(I).LT.ENERGY(ISTATE)) ISTATE=I
       END DO
       DO I=1,NSTATE
        RELAU = ENERGY(I)-ENERGY(ISTATE)
        RELEV = RELAU * CONV_AU_TO_EV_
        RELCM = RELAU * CONV_AU_TO_CM1_
        RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
        WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &   MSTATE(I), RELAU, RELEV, RELCM, RELKJ
       END DO
       WRITE(6,*)
      END IF

      IF(NLYROOT.NE.0) IFMSCOUP=.FALSE.
      IF(IFMSCOUP) THEN
        Call StatusLine('CASPT2:','Effective Hamiltonian')
        CALL MLTCTL(HEFF,UEFF,U0)

        IF(IPRGLB.GE.VERBOSE.AND.(NLYROOT.EQ.0)) THEN
         WRITE(6,*)' Relative (X)MS-CASPT2 energies:'
         WRITE(6,'(1X,A4,4X,A12,1X,A10,1X,A10,1X,A10)')
     &     'Root', '(a.u.)', '(eV)', '(cm^-1)', '(kJ/mol)'
         DO I=1,NSTATE
          RELAU = ENERGY(I)-ENERGY(1)
          RELEV = RELAU * CONV_AU_TO_EV_
          RELCM = RELAU * CONV_AU_TO_CM1_
          RELKJ = RELAU * CONV_AU_TO_KJ_PER_MOLE_
          WRITE(6,'(1X,I4,4X,F12.8,1X,F10.2,1X,F10.1,1X,F10.2)')
     &     I, RELAU, RELEV, RELCM, RELKJ
         END DO
         WRITE(6,*)
        END IF
      END IF

* Back-transform the effective Hamiltonian and the transformation matrix
* to the basis of original CASSCF states
      CALL Backtransform(Heff,Ueff,U0)

* create a JobMix file
* (note that when using HDF5 for the PT2 wavefunction, IFMIX is false)
      CALL CREIPH_CASPT2(Heff,Ueff,U0)

* Store the PT2 energy and effective Hamiltonian on the wavefunction file
      CALL PT2WFN_ESTORE(HEFF)

* Store rotated states if XMUL + NOMUL
      IF (IFXMS.AND.(.NOT.IFMSCOUP)) CALL PT2WFN_DATA

* store information on runfile for geometry optimizations
      Call Put_iScalar('NumGradRoot',iRlxRoot)
      Call Store_Energies(NSTATE,ENERGY,iRlxRoot)

      CALL MMA_DEALLOCATE(UEFF)
      CALL MMA_DEALLOCATE(U0)
      End If  !Skipping MultiState calculation when IFNOPT2=true
9000  CONTINUE

C Free resources, close files
      CALL PT2CLS
C     If (nRoots.ne.lRoots) Call ModDip

      CALL MMA_DEALLOCATE(HEFF)
      CALL MMA_DEALLOCATE(H0)

C     PRINT I/O AND SUBROUTINE CALL STATISTICS
      IF ( IPRGLB.GE.USUAL ) THEN
        CALL FASTIO('STATUS')
        CALL QSTAT(' ')
      END IF

      Call StatusLine('CASPT2:','Finished.')
      CALL QEXIT('CASPT2')
      RETURN
      END
      subroutine getfdiag(iorb,eig,nbas)
      implicit real*8 (a-h,o-z)
      dimension tmp(nbas)
      Call Get_dArray('RASSCF OrbE',tmp,nbas)
C     if (iorb.eq.1) eig = -20.4907280456d+00
C     if (iorb.eq.2) eig = -11.1340416978D+00
      eig=tmp(iorb)
      end
