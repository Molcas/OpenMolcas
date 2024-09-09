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
      SUBROUTINE FOCK_m(F,FI,FP,D,P,Q,FINT)
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
      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use rasscf_global, only: E2act, ECAS, HALFQ1,
     &                         NTOT3, ISTORP,
     &                         ISTORD, ITRI, CBLB, iBLB,
     &                         jBLB
      IMPLICIT None

      REAL*8 F(*),FP(*),D(*),P(*),Q(*),FINT(*),FI(*)
      integer ISTSQ(8),ISTAV(8)
      real*8 ECAS0

#include "rasdim.fh"
#include "general.fh"
      Character(LEN=16), Parameter:: ROUTINE='FOCK    '
      Integer iPrLev
      REAL*8  CSX, QNTM
      Integer ipFMCSCF, ISTBM, ISTD, ISTFCK, ISTFP, ISTP, ISTZ,
     &        ISYM, IX, IX1, JSTF, N1, N2, NAO, NEO, NI, NIA, NIO,
     &        NM, NO, NO2, NOR, NP, NT, NTM, NTT,
     &        NTV, NUVX, NV, NVI, NVM

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
C *****************************************

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
      E2act=0.0d0
      HALFQ1=0.0D0
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
       CALL FZERO(F(ISTFCK+1),NO**2)
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


c
c          first compute the Q-matrix (equation (19))
c
c          Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c
c          P is packed in xy and pre-multiplied by 2
c                            and reordered
c
        CALL DGEMM_('N','N',
     &              NO,NAO,NUVX,
     &              1.0d0,FINT(JSTF),NO,
     &              P(ISTP),NUVX,
     &              0.0d0,Q,NO)


c
c       active-active interaction term in the RASSCF energy
c
        ECAS0=ECAS
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
!          ECAS=ECAS+0.5D0*Q(NTT)
         HALFQ1=HALFQ1+0.5D0*Q(NTT)
         E2act=E2act+0.5D0*Q(NTT)
        END DO


!       Fock matrix
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

! Unknown below

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

! End of long loop over symmetry
      END DO

      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A)')' MCSCF Fock-matrix in MO-basis'
        ipFMCSCF=1
        Do iSym=1,nSym
           nOr=nOrb(iSym)
           Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
           ipFMCSCF=ipFMCSCF+nOr*nOr
        End Do
      End If

      END SUBROUTINE FOCK_m

      SUBROUTINE FOCK_update(F,FI,FP,D,P,Q,FINT,CMO)
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
C          ********** IBM-3090 MOLCASs Release: 90 02 22 **********
C
      use printlevel, only: debug
      use mspdft, only: iIntS
      use mcpdft_output, only: lf, iPrLoc
      use mspdftgrad,only:FxyMS
      use mcpdft_input, only: mcpdft_options
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: E2act, nTot3, nTot4, ISTORP, ISTORD,
     &                         iTri, CBLB, IBLB, JBLB

      IMPLICIT None
      REAL*8 FI(*),FP(*),D(*),P(*),Q(*),FINT(*),F(*),CMO(*)
      integer ISTSQ(8),ISTAV(8)

      Real*8, Allocatable:: TF(:)

#include "rasdim.fh"
#include "general.fh"
      Character(LEN=16), Parameter:: ROUTINE='FOCK    '
      Integer iPrLev
      Real*8 CSX, E2eP, QNTM
      Integer i, ipFMCSCF, ISTBM, ISTD, ISTFCK, ISTFP, ISTP, ISTZ,
     &        iSym, IX1, JSTF, N1, N2, NAO, NEO, NI, NIO, NM, NO, NO2,
     &        NOR, NP, NT, NTM, NTT, NTV, NUVX, NV, NVI, NVM

C
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      Call mma_allocate(TF,NTOT4,Label='TF')
      TF(:)=0.0D0


      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do
! *****************************************


!     add FI to FA to obtain FP
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)
!     LOOP OVER ALL SYMMETRY BLOCKS

      ISTFCK=0
      ISTFP=0
      ISTD=0
      ISTBM=0
      IX1=0
      ISTZ=0
      E2act=0.0d0

! A long loop over symmetry
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NEO=NSSH(ISYM)
       NO=NORB(ISYM)
       NO2=(NO**2+NO)/2
       CSX=0.0D0
       N1=0
       N2=0
       IF(NO.EQ.0) GO TO 90
       CALL FZERO(TF(ISTFCK+1),NO**2)

!    First index in F is inactive

       IF(NIO.NE.0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          TF(ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
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


c
c          first compute the Q-matrix (equation (19))
c
c          Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c
c          P is packed in xy and pre-multiplied by 2
c                            and reordered

        CALL DGEMM_('N','N',
     &              NO,NAO,NUVX,
     &              1.0d0,FINT(JSTF),NO,
     &              P(ISTP),NUVX,
     &              0.0d0,Q,NO)


!Now Q should contain the additional 2-electron part of the fock matrix
!for mcpdft, for the active region, at least.

!We should also have contributions from terms like FI and FA, too.
!FA takes care of the 1-RDM/2e- integral terms?
!FI takes care of the one-body hamiltonian and the occ/occ and occ/act
!contributions.

        E2eP=0d0
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
         E2eP=E2eP+0.5D0*Q(NTT)
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
          TF(ISTFCK+NO*(NM-1)+NT+NIO)=QNTM
         END DO
        END DO
       ENDIF


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
        write(6,*) TF(i)
      end do
      call xflush(6)
      end if
      Call DAXPY_(NTOT4,1.0d0,TF,1,F,1)
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

!For MCLR
      IF(mcpdft_options%grad .and. mcpdft_options%mspdft) THEN
       CALL DCopy_(nTot4,F,1,FxyMS(:,iIntS),1)
      ELSE
       Call put_dArray('Fock_PDFT',F,ntot4)
      END IF

      CALL FOCKOC(Q,F,CMO)

      Call mma_deallocate(TF)

      END SUBROUTINE FOCK_update
