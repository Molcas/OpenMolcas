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
      SUBROUTINE FOCK_m(F,BM,FI,FP,D,P,Q,FINT,IFINAL,CMO)
C
C     RASSCF program version IBM-3090: SX section
c
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
C
C *** Cholesky section ********************
c      Call DecideOnCholesky(DoCholesky)

      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do
C *****************************************


      ipFint = ip_Dummy
      ipP2reo= ip_Dummy
      If(KSDFT(1:3).ne.'SCF'.and.
     &         DFTFOCK(1:4).eq.'DIFF'.and.
     &         nac.ne.0) Then
        Call GetMem('TmpPUVX','Allo','Real',ipFint,nFint)
c       Call Get_Temp('TmpPUVX ',Work(ipFint),nFint)
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
c
c     add FI to FA to obtain FP
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)
C     LOOP OVER ALL SYMMETRY BLOCKS
C
      ISTFCK=0
      ISTFP=0
      ISTD=0
      ISTBM=0
      IX1=0
      ISTZ=0
      ioffQmat=0
      E2act=0.0d0
C
* A long loop over symmetry
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
       IF(NO.EQ.0) GO TO 90
       CALL VCLR(F(ISTFCK+1),1,NO**2)
c
c      first index in F is inactive
c
       IF(NIO.NE.0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          F(ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
         END DO
        END DO
       ENDIF
c
c      first index in F active
c
       IF(NAO.NE.0) THEN

        ISTP=ISTORP(ISYM)+1
        JSTF=ISTORD(ISYM)+1
        NUVX=(ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

        If (.not.DoCholesky .or. ALGO.eq.1) Then
c
c          first compute the Q-matrix (equation (19))
c
c          Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c
c          P is packed in xy and pre-multiplied by 2
c                            and reordered
c
c          write(6,*) 'PUVX integrals in FOCK'
c         call wrtmat(FINT(JSTF),1,nFInt,1,nFInt)
c         write(6,*) 'two-elec density mat OR DMAT*DMAT in FOCK'
c         call wrtmat(P(ISTP),1,nFint,1,nFint)
           CALL DGEMM_('N','N',
     &                 NO,NAO,NUVX,
     &                 1.0d0,FINT(JSTF),NO,
     &                 P(ISTP),NUVX,
     &                 0.0d0,Q,NO)

        ElseIf (ALGO.eq.2) Then

c --- the Q-matrix has been already computed as Q(a,v)
c --- where a is an AO index and v is an active index
c ---
c --- Transform the 1st index to MOs (one symmetry at the time)
c ---
c ---     Q(m,v) = C(a,m) * Q(a,v)

          ipQS = ipQmat + ISTAV(iSym)
          ipMOs= 1 + ISTSQ(iSym) + nBas(iSym)*nFro(iSym)

          CALL DGEMM_('T','N',nOrb(iSym),nAsh(iSym),nBas(iSym),
     &               1.0d0,CMO(ipMOs),nBas(iSym),
     &               Work(ipQS),nBas(iSym),
     &               0.0d0,Q(1),nOrb(iSym))

          write(6,*) 'transforming the Q-matrix'

        Else

          Write(LF,*)'FOCK: illegal Cholesky parameter ALGO= ',ALGO
          call qtrace()
          call abend()

        EndIf

CGLM        call recprt('Q-mat',' ',Q(1),NO,NAO)

c
c       active-active interaction term in the RASSCF energy
c
        ECAS0=ECAS
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
         ECAS=ECAS+0.5D0*Q(NTT)
         HALFQ1=HALFQ1+0.5D0*Q(NTT)
         E2act=E2act+0.5D0*Q(NTT)
        END DO
      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'Two-electron contribution (Q term):', ECAS-ECAS0
      END IF
C        write(6,*) 'ECAS aft adding Q in FOCK :', ECAS
*
        If(ipFint.ne.ip_Dummy) Then
          Call GetMem('TmpQ','Allo','Real',ipQ,NAO*NO)
          If (ipP2reo.ne.ip_Dummy) Then
             CALL DGEMM_('N','N',
     &                   NO,NAO,NUVX,
     &                   1.0d0,Work(ipFint+JSTF-1),NO,
     &                   Work(ipP2reo+ISTP-1),NUVX,
     &                   0.0d0,Work(ipQ),NO)
          Else
             CALL DGEMM_('N','N',
     &                   NO,NAO,NUVX,
     &                   1.0d0,Work(ipFint+JSTF-1),NO,
     &                   P(ISTP),NUVX,
     &                   0.0d0,Work(ipQ),NO)
          End If
          Call DaXpY_(NAO*NO,1.0d0,Work(ipQ),1,Q,1)
*
          DO NT=1,NAO
            NTT=(NT-1)*NO+NIO+NT
            HALFQ=HALFQ+0.5D0*Work(ipQ+NTT-1)
          END DO
          Call GetMem('TmpQ','Free','Real',ipQ,NAO*NO)
        End If

c
c       Fock matrix
c
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
c
c       The Brillouin matrix BM(pq)=F(qp)-F(pq)
c
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
c
* End of long loop over symmetry
      END DO
c
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
C
c     Calculate Fock matrix for occupied orbitals.
C
      If (iFinal.eq.1) CALL FOCKOC_m(Q,F,CMO)
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

      SUBROUTINE FOCK_update(F,BM,FI,FP,D,P,Q,FINT,IFINAL,CMO)
!This subroutine is supposed to add the dft portions of the mcpdft fock
!matrix to the Fock matrix pieces that have already been built for the
!CASSCF portion.

C
C     RASSCF program version IBM-3090: SX section
c
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
      integer ISTSQ(8),ISTAV(8),iTF

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

      Call Unused_real_array(BM)
      Call Unused_integer(ifinal)
      Call GetMem('fockt','ALLO','REAL',iTF,NTOT4)
      Call dcopy_(ntot4,0d0,0,Work(iTF),1)
C
C *** Cholesky section ********************
c      Call DecideOnCholesky(DoCholesky)

      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do
C *****************************************

!      Call GetMem('ONTOPT','ALLO','Real',iTEOTP,NFINT)
!      Call GetMem('ONTOPO','ALLO','Real',iOEOTP,NTOT1)
      !Read in the one- and two- electron potentials.
!      Call Get_dArray('ONTOPT',work(iTEOTP),NFINT)
!      Call Get_dArray('ONTOPO',work(iOEOTP),NTOT1)


!I think the best way forward is to construct FI, FA (MO basis) using
!the potentials (v_pqrs and V_pq) instead of the integrals.  I think we
!want to use the full 2-body density matrix and ExFac = 1.  If we have
!FI, FA, and Q, then we can use the prescription from fock.f to
!construct the Focc term (the part of the fock matrix that we want).



!******************************************************************
!
! Build the FA terms using the potentials v.
!
!******************************************************************

!      ExFac_tmp = 1.0d0
!      Call Upd_FA_m(Work(iTEOTP),FP,D,ExFac_tmp)
!Check - does this regenerate FA if the regular integrals are passed?
!FP should contain the Fock matrix contribution that we want.


!******************************************************************
!
! Build the FI terms using the potentials v and V.
!
!******************************************************************


!      iOff1 = 0
!      Do ISYM=1,NSYM
!        Do iOrb=1,norb(iSym)
!          do jOrb=1,iOrb
!should we follow the guide of ftwo.f?
!for starters, I don't seem to have all the necessary 2-body potentials,
!right?

!          end do
!        end do
!      end do


c     add FI to FA to obtain FP
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)
C     LOOP OVER ALL SYMMETRY BLOCKS

      ISTFCK=0
      ISTFP=0
      ISTD=0
      ISTBM=0
      IX1=0
      ISTZ=0
      ioffQmat=0
      E2act=0.0d0
C
* A long loop over symmetry
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
       IF(NO.EQ.0) GO TO 90
       CALL VCLR(Work(iTF-1+ISTFCK+1),1,NO**2)

!    First index in F is inactive

       IF(NIO.NE.0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          Work(iTF-1+ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
         END DO
        END DO
       ENDIF
c
c      first index in F active
c
       IF(NAO.NE.0) THEN

        ISTP=ISTORP(ISYM)+1
        JSTF=ISTORD(ISYM)+1
        NUVX=(ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

        If (.not.DoCholesky .or. ALGO.eq.1) Then
c
c          first compute the Q-matrix (equation (19))
c
c          Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c
c          P is packed in xy and pre-multiplied by 2
c                            and reordered
c
c          write(6,*) 'PUVX integrals in FOCK'
c         call wrtmat(FINT(JSTF),1,nFInt,1,nFInt)
c         write(6,*) 'two-elec density mat OR DMAT*DMAT in FOCK'
c         call wrtmat(P(ISTP),1,nFint,1,nFint)
           CALL DGEMM_('N','N',
     &                 NO,NAO,NUVX,
     &                 1.0d0,FINT(JSTF),NO,
     &                 P(ISTP),NUVX,
     &                 0.0d0,Q,NO)


!Now Q should contain the additional 2-electron part of the fock matrix
!for mcpdft, for the active region, at least.

!We should also have contributions from terms like FI and FA, too.
!FA takes care of the 1-RDM/2e- integral terms?
!FI takes care of the one-body hamiltonian and the occ/occ and occ/act
!contributions.

!        write(*,*) 'q-matrix'
!        do i=1,nO*nAO
!          write(*,*) Q(i)
!        end do

        Else

          Write(LF,*)'FOCK: illegal Cholesky parameter ALGO= ',ALGO
          call qtrace()
          call abend()

        EndIf

        E2eP=0d0
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
         E2eP=E2eP+0.5D0*Q(NTT)
         !ECAS=ECAS+Q(NTT)
        END DO
c
c       Fock matrix
c
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
          Work(iTF-1+ISTFCK+NO*(NM-1)+NT+NIO)=QNTM
         END DO
        END DO
       ENDIF

CGLM        call recprt('Q-mat',' ',Q(1),NO,NAO)

c
c       active-active interaction term in the RASSCF energy
c
c
* End of long loop over symmetry
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


!Now, add all components to the original Fock matrix
!      Call DAXPY_()
!      Call DAXPY_()
c
C
c     Calculate Fock matrix for occupied orbitals.
C

      If ( iPrLev.ge.DEBUG ) then
      write(6,*) 'old fock terms:'
      do i=1,Ntot4
        write(6,*) F(i)
      end do
      write(6,*) 'new fock terms to add:'
      do i=1,Ntot4
        write(6,*) Work(itF-1+i)
      end do
      call xflush(6)
      end if
      Call DAXPY_(NTOT4,1.0d0,Work(iTF),1,F,1)
!      write(*,*) 'added new fock terms to old fock matrix'
!I am going to add the Fock matrix temporarily to the Runfile.  I don't
!want to construct it again in MCLR in the case of gradients.
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A)')' MCSCF Fock-matrix in MO-basis'
        ipFMCSCF=1
        Do iSym=1,nSym
           nOr=nOrb(iSym)
           Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
           ipFMCSCF=ipFMCSCF+nOr*nOr
        End Do
      End If

!What happens if we divide by two?
!      Call Dscal_(ntot4,0.5d0,F,1)

!For MCLR
      Call put_dArray('Fock_PDFT',F,ntot4)

      call xflush(6)
      CALL FOCKOC_m(Q,F,CMO)
C
      Call GetMem('fockt','Free','REAL',iTF,NTOT4)
C
      RETURN
      END
