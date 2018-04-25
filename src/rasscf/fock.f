************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE FOCK(F,BM,FI,FP,D,P,Q,FINT,IFINAL,CMO)
C
C     RASSCF program version IBM-3090: SX section
c
******************************************************************************
*   BM    : output
*   F     : used only in this routine
*   Q     : used only in this routine
*   FP    : in Input it is the FA matrix in MO. in Output it is FI + FA in MO.
*   FI    : input. FI in MO basis
*   D     : input. Active 1RDM in MO basis
*   P     : input. Active 2RDM in MO basis
*   FINT  : input. Active two-electron integrals in MO basis.
*   IFINAL: input
*   CMO   : input
******************************************************************************
c     Calculation of the MCSCF fock matrix F(eq.(7) in I.J.Q.C.S14,175)
c     FP is the matrix FI+FA (FP is FA at entrance)
c     F is stored as a symmetry blocked square matrix, by columns.
c     Note that F contains all elements, also the zero elements
c     occurring when the first index is secondary.
c     F is used to construct the Brillouin elements and the first row
c     of the super-CI Hamiltonian, while FP is used as the effective
c     one-electron operator in the construction of the super-CI
c     interaction matrix.
c
C          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FI(*),FP(*),D(*),P(*),Q(*),FINT(*),F(*),BM(*),CMO(*)
      integer ISTSQ(8),ISTAV(8)
      real*8 ECAS0

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='FOCK    ')
#include "WrkSpc.fh"
      Logical DoActive,DoQmat,DoCholesky
      Integer ALGO

      COMMON /CHOTODO /DoActive,DoQmat,ipQmat
      COMMON /CHLCAS  /DoCholesky,ALGO
C
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do

      ipFint = ip_Dummy
      ipP2reo= ip_Dummy
      If(KSDFT(1:3).ne.'SCF'.and.
     &         DFTFOCK(1:4).eq.'DIFF'.and.
     &         nac.ne.0) Then
        Call GetMem('TmpPUVX','Allo','Real',ipFint,nFint)
        Call Get_dArray('DFT_TwoEl',Work(ipFint),nFint)
        HALFQ=0.0d0
        HALFQ1=0.0d0
        If(Exfac.ne.1.0d0) Then
          Call Get_Temp('nP2reo  ',P2reo,1)
          nP2reo=Int(P2reo)
          CALL GETMEM('P2_reo','ALLO','REAL',ipP2reo,nP2reo)
          Call Get_Temp('P2_reo  ',Work(ipP2reo),nP2reo)
        End If
      End If

*********************************************************************************
* add FI to FA to obtain FP
*********************************************************************************
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)

*********************************************************************************
* Loop over symmetry blocks
*********************************************************************************
      ISTFCK=0
      ISTFP=0
      ISTD=0
      ISTBM=0
      IX1=0
      ISTZ=0
      ioffQmat=0
      E2act=0.0d0

      DO ISYM=1,NSYM
       IX=IX1+NFRO(ISYM)
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NEO=NSSH(ISYM)
       NIA=NIO+NAO
       NO=NORB(ISYM)
       NO2=(NO**2+NO)/2
       CSX=0.0D0
       N1=0
       N2=0

* If no orbitals in this symmetry block go to next symm
       IF(NO.EQ.0) GO TO 90

*************************
* clear F_gen matrix
*************************
       CALL FZERO(F(ISTFCK+1),NO**2)

*********************************************************************************
* first index in F is inactive: F_gen is twice FP=FI+FA (Eq.10.8.27 MEST)
*********************************************************************************
       IF(NIO.NE.0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          F(ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
         END DO
        END DO
       ENDIF

*********************************************************************************
* first index in F active
*********************************************************************************
       IF(NAO.NE.0) THEN
        ISTP=ISTORP(ISYM)+1
        JSTF=ISTORD(ISYM)+1
        NUVX=(ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

        If (.not.DoCholesky .or. ALGO.eq.1) Then
*********************************************************************************
* Compute the Q-matrix (Eq. (19) in IJQC S14 175 1980 or Eq. 10.8.31 MEST)
* Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
* P is packed in xy and pre-multiplied by 2 and reordered
*********************************************************************************
           IF(IPRLEV.ge.DEBUG) THEN
             write(6,*) 'PUVX integrals in FOCK'
             call wrtmat(FINT(JSTF),NO,NUVX,NO,NUVX)
             write(6,*) 'two-elec density mat OR DMAT*DMAT in FOCK'
             call wrtmat(P(ISTP),NAO,NUVX,NAO,NUVX)
           END IF

           CALL DGEMM_('N','N',
     &                 NO,NAO,NUVX,
     &                 1.0d0,FINT(JSTF),NO,
     &                 P(ISTP),NUVX,
     &                 0.0d0,Q,NO)

        ElseIf (ALGO.eq.2) Then
*********************************************************************************
c --- the Q-matrix has been already computed as Q(a,v)
c --- where a is an AO index and v is an active index
c --- Transform the 1st index to MOs (one symmetry at the time)
c --- Q(m,v) = C(a,m) * Q(a,v)
*********************************************************************************
          ipQS = ipQmat + ISTAV(iSym)
          ipMOs= 1 + ISTSQ(iSym) + nBas(iSym)*nFro(iSym)

          CALL DGEMM_('T','N',nOrb(iSym),nAsh(iSym),nBas(iSym),
     &               1.0d0,CMO(ipMOs),nBas(iSym),
     &               Work(ipQS),nBas(iSym),
     &               0.0d0,Q(1),nOrb(iSym))
        Else
          Write(LF,*)'FOCK: illegal Cholesky parameter ALGO= ',ALGO
          call qtrace()
          call abend()
        EndIf

        IF(IPRLEV.ge.DEBUG) THEN
          call recprt('Q-mat in fock.f',' ',Q(1),NO,NAO)
        END IF

*********************************************************************************
* active-active interaction energy term in the RASSCF energy: trace of Q
*********************************************************************************
        ECAS0=ECAS
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
         ECAS=ECAS+0.5D0*Q(NTT)
         HALFQ1=HALFQ1+0.5D0*Q(NTT)
         E2act=E2act+0.5D0*Q(NTT)
        END DO
        IF(IPRLEV.ge.DEBUG) THEN
          write(6,*) 'Two-electron contribution (Q term):', ECAS-ECAS0
          write(6,*) 'ECAS aft adding Q in fock.f :', ECAS
        END IF

*********************************************************************************
* Is this relevant?
          If(ipFint.ne.ip_Dummy) Then
            Call GetMem('TmpQ','Allo','Real',ipQ,NAO*NO)
            If (ipP2reo.ne.ip_Dummy) Then
               CALL DGEMM_('N','N',
     &                     NO,NAO,NUVX,
     &                     1.0d0,Work(ipFint+JSTF-1),NO,
     &                     Work(ipP2reo+ISTP-1),NUVX,
     &                     0.0d0,Work(ipQ),NO)
            Else
               CALL DGEMM_('N','N',
     &                     NO,NAO,NUVX,
     &                     1.0d0,Work(ipFint+JSTF-1),NO,
     &                     P(ISTP),NUVX,
     &                     0.0d0,Work(ipQ),NO)
            End If
            Call DaXpY_(NAO*NO,1.0d0,Work(ipQ),1,Q,1)
*
            DO NT=1,NAO
              NTT=(NT-1)*NO+NIO+NT
              HALFQ=HALFQ+0.5D0*Work(ipQ+NTT-1)
            END DO
            Call GetMem('TmpQ','Free','Real',ipQ,NAO*NO)
          End If

*********************************************************************************
* Continue... Fock matrix for first index active
*********************************************************************************
          NTM=0
          DO NT=1,NAO
           DO NM=1,NO
            NTM=NTM+1
            QNTM=Q(NTM)
            DO NV=1,NAO
             NVI=NV+NIO
             NTV=ITRI(MAX(NT,NV))+MIN(NT,NV)+ISTD
             NVM=ITRI(MAX(NVI,NM))+MIN(NVI,NM)+ISTFP
             QNTM=QNTM+D(NTV)*FI(NVM)
            END DO
            F(ISTFCK+NO*(NM-1)+NT+NIO)=QNTM
           END DO
          END DO

       ENDIF
*********************************************************************************
*      ^ End if related to first index in F active
*********************************************************************************

*********************************************************************************
*       The Brillouin matrix BM(pq)=F(qp)-F(pq)
*********************************************************************************
       NPQ=ISTBM
       DO NP=NIO+1,NO
        DO NQ=1,NIA
         NPQ=NPQ+1
         BM(NPQ)=F(ISTFCK+NO*(NP-1)+NQ)-F(ISTFCK+NO*(NQ-1)+NP)
c
c        Set zeroes to BM elements corresponding to rotations not allowed
c        as controlled by the array IZROT.
c
         IF(NP.LE.NIA.AND.NQ.GT.NIO) THEN
          NT=NP-NIO
          NU=NQ-NIO
          IF(NT.LE.NU) THEN
           BM(NPQ)=0.0D0
          ELSE
           NTU=ISTZ+ITRI(NT-1)+NU
           IF(IZROT(NTU).NE.0) THEN
            BM(NPQ)=0.0D0
C            Write(LF,*)'FOCK: IZROT=1 so BLB=0 for NP,NQ=',NP,NQ
           END IF
           IF(IXSYM(IX+NP).NE.IXSYM(IX+NQ)) THEN
            BM(NPQ)=0.0D0
C            Write(LF,*)'FOCK: IXSYM=1 so BLB=0 for NP,NQ=',NP,NQ
           END IF
          ENDIF
         ENDIF
c
c        check for largest Brillouin matrix element
c
         IF(ABS(BM(NPQ)).LT.ABS(CSX)) GO TO 20
         IF(IXSYM(IX+NP).NE.IXSYM(IX+NQ)) GO TO 20
         CSX=BM(NPQ)
         N1=NQ+NFRO(ISYM)
         N2=NP+NFRO(ISYM)
 20      CONTINUE
        END DO
       END DO
c
90     CONTINUE
       ISTFCK=ISTFCK+NO**2
       ISTFP=ISTFP+NO2
       ISTD=ISTD+(NAO**2+NAO)/2
       ISTBM=ISTBM+(NIO+NAO)*(NAO+NEO)
       IX1=IX1+NBAS(ISYM)
       ISTZ=ISTZ+(NAO**2-NAO)/2
       CBLB(ISYM)=CSX
       IBLB(ISYM)=N1
       JBLB(ISYM)=N2
      END DO
*********************************************************************************
*     ^ End of long loop over symmetry
*********************************************************************************

      If ( iPrLev.ge.DEBUG ) then
        CASDFT_En=0.0d0
        If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM')
     &   Call Get_dScalar('CASDFT energy',CASDFT_En)
        Write(LF,'(A,2F22.16)') ' RASSCF energy: ',
     &                  ECAS+CASDFT_En,VIA_DFT
      End If
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A)')' MCSCF Fock-matrix in MO-basis'
        ipFMCSCF=1
        Do iSym=1,nSym
           nOr=nOrb(iSym)
           Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
           ipFMCSCF=ipFMCSCF+nOr*nOr
        End Do
      End If
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A)')' Brillouin matrix'
        ipBM=1
        Do iSym=1,nSym
           nIs=nIsh(iSym)
           nAs=nAsh(iSym)
           nSs=nSsh(iSym)
           Call RecPrt(' ',' ',BM(ipBM),(nAs+nIs),(nAs+nSs))
           ipBM=ipBM+(nAs+nIs)*(nAs+nSs)
        End Do
      End If
C
C     Maximum BLB matrix element all symmetries
C
      CBLBM=0.0D0
      ISYMBB=0
      DO ISYM=1,NSYM
       IF(ABS(CBLB(ISYM)).LT.ABS(CBLBM)) GO TO 150
       CBLBM=CBLB(ISYM)
       IBLBM=IBLB(ISYM)
       JBLBM=JBLB(ISYM)
       ISYMBB=ISYM
 150   CONTINUE
      END DO

*********************************************************************************
*     Calculate Fock matrix for occupied orbitals.
*********************************************************************************
      If (iFinal.eq.1) CALL FOCKOC(Q,F,CMO)
C
      If(ipFint.ne.ip_Dummy) Then
        Call GetMem('TmpPUVX','Free','Real',ipFint,nFint)
      End If

      If(ipP2reo.ne.ip_Dummy) Then
        CALL GETMEM('P2_reo','FREE','REAL',ipP2reo,nP2reo)
      End If

      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' >>> Exit Fock <<< '
         Write(LF,*)
      End If
C
      RETURN
      END
