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
      SUBROUTINE TSHinit(Energy)
      use rasdef, only: NRAS, NRASEL, NRSPRT, NRS1, NRS1T, NRS2, NRS3
      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: PART, JBNUM, LROOT
      use gugx, only: SGStruct, CIStruct, EXStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, LSYM1, LSYM2, IRREP, MLTPLT,
     &                 NACTE, NELE3, NHOLE1, RASTYP
      use cntrl, only: ISTATE1, nCI1, ISTATE2, nCI2, ChkHop
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NDEL,NFRO,NISH,NSSH

      IMPLICIT None
      Real*8 :: Energy(nState)

      Type (SGStruct), Target :: SGS(2)
      Type (CIStruct) :: CIS(2)
      Type (EXStruct) :: EXS(2)
      INTEGER      I,JOB1,JOB2,iRlxRoot
      CHARACTER(LEN=8) WFTYP1,WFTYP2
      LOGICAL   LOWROOT, UPROOT
      Real*8, Allocatable:: CI1(:), CI2(:)
      Integer NACTE1, MPLET1, NHOL11, NELE31, NACTE2, MPLET2, NHOL12,
     &        NELE32
      REAL*8 EDIFF, EThr

      Interface
      Subroutine SGInit(nSym,nActEl,iSpin,SGS,CIS)
      use gugx, only: SGStruct, CIStruct
      IMPLICIT None
      Integer nSym, nActEl, iSpin
      Type (SGStruct), Target :: SGS
      Type (CIStruct) :: CIS
      End Subroutine SGInit
      End Interface

*
C
C Print a banner
C
      IF (IPGLOB.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(6X,A)') repeat('*',100)
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A,36X,A,37X,A)')
     &       '*',' Surface hopping section ','*'
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A)') repeat('*',100)
        WRITE(6,*)
        WRITE(6,*)

        WRITE(6,'(6X,A)')'Surface hopping section'
        WRITE(6,'(6X,A)')'-----------------------'
      END IF
C
C Get the current state and print it's energy
C
      CALL Get_iScalar('Relax CASSCF root',iRlxRoot)
      IF (IPGLOB.GE.2) THEN
        WRITE(6,'(6X,A,I14)')'The current state is:',iRlxRoot
        WRITE(6,'(6X,A,6X,ES15.6,A,/)')'Its energy is:',
     &                                 ENERGY(iRlxRoot),' a.u.'
      END IF
C
C Get wave function parameters for current state
C
      ISTATE1=iRlxRoot
      JOB1=JBNUM(ISTATE1)
      NACTE1=NACTE(JOB1)
      MPLET1=MLTPLT(JOB1)
      LSYM1=IRREP(JOB1)
      NHOL11=NHOLE1(JOB1)
      NELE31=NELE3(JOB1)
      WFTYP1=RASTYP(JOB1)
      SGS(1)%IFRAS=1
      SGS(2)%IFRAS=1
C
C Set the variables for the wave function of the current state
C
      IF(WFTYP1.EQ.'GENERAL ') THEN
          NRSPRT=3
          DO I=1,8
             NRAS(I,1)=NRS1(I)
             NRAS(I,2)=NRS2(I)
             NRAS(I,3)=NRS3(I)
          END DO
          NRASEL(1)=2*NRS1T-NHOL11
          NRASEL(2)=NACTE1-NELE31
          NRASEL(3)=NACTE1
          CALL SGINIT(NSYM,NACTE1,MPLET1,SGS(1),CIS(1))
          IF(IPGLOB.GT.4) THEN
             WRITE(6,*)'Split-graph structure for JOB1=',JOB1
             CALL SGPRINT(SGS(1))
          END IF
          CALL CXINIT(SGS(1),CIS(1),EXS(1))
C CI sizes, as function of symmetry, are now known.
          NCI1=CIS(1)%NCSF(LSYM1)
      ELSE
* The only other cases are HISPIN, CLOSED or EMPTY.
* NOTE: The HISPIN case is suspected to be buggy. Not used now.
          NCI1=1
      END IF
      CALL mma_allocate(CI1,NCI1,Label='CI1')

C
C Check for the possibility to hop to a lower root
C
      IF ((ISTATE1-1).GE.1) THEN
         LOWROOT=.TRUE.
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I2)') "There is a lower root, which is: "
     &     ,LROOT(ISTATE1-1)
         END IF
      ELSE
         LOWROOT=.FALSE.
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A)') "There is no lower root"
         END IF
      END IF
C
C Check for the possibility to hop to an upper root
C
      IF ((ISTATE1+1).LE.NSTATE) THEN
         UPROOT=.TRUE.
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I2,/)') "There is an upper root, which is: "
     &     ,LROOT(ISTATE1+1)
         END IF
      ELSE
         UPROOT=.FALSE.
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,/)') "There is no upper root"
         END IF
      END IF
C
C Check surface hopping to a root lower than the current one
C
      IF (LOWROOT) THEN
         ISTATE2=ISTATE1-1
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I3)') 'The lower state is:',ISTATE2
           WRITE(6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',
     &          ENERGY(ISTATE2),' a.u.'
           WRITE(6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',
     &          ENERGY(iRlxRoot)-ENERGY(ISTATE2),' a.u.'
         END IF
*---------------------------------------------------------------------*
* Adapted from RASSI subroutines GTMCTL and READCI                    *
*                                                                     *
C Get wave function parameters for ISTATE2
         JOB2=JBNUM(ISTATE2)
         NACTE2=NACTE(JOB2)
         MPLET2=MLTPLT(JOB2)
         LSYM2=IRREP(JOB2)
         NHOL12=NHOLE1(JOB2)
         NELE32=NELE3(JOB2)
         WFTYP2=RASTYP(JOB2)
         Call NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
         IF(IPGLOB.GE.4) CALL PRPRTTAB(PART)
C For the second wave function
         IF(WFTYP2.EQ.'GENERAL ') THEN
            NRSPRT=3
            DO I=1,8
               NRAS(I,1)=NRS1(I)
               NRAS(I,2)=NRS2(I)
               NRAS(I,3)=NRS3(I)
            END DO
            NRASEL(1)=2*NRS1T-NHOL12
            NRASEL(2)=NACTE2-NELE32
            NRASEL(3)=NACTE2
            CALL SGINIT(NSYM,NACTE2,MPLET2,SGS(2),CIS(2))
            IF(IPGLOB.GT.4) THEN
               WRITE(6,*)'Split-graph structure for JOB2=',JOB2
               CALL SGPRINT(SGS(2))
            END IF
            CALL CXINIT(SGS(2),CIS(2),EXS(2))
C     CI sizes, as function of symmetry, are now known.
            NCI2=CIS(2)%NCSF(LSYM2)
         ELSE
C     Presently, the only other cases are HISPIN, CLOSED or EMPTY.
            NCI2=1
         END IF
         CALL mma_allocate(CI2,NCI2,Label='CI2')
C     Check if the Energy gap is smaller than the threshold.
         Ediff=ABS(ENERGY(ISTATE2)-ENERGY(ISTATE1))
         Ethr=3.0D-2
         IF (Ediff.le.Ethr) THEN
            ChkHop=.TRUE.
         ELSE
            ChkHop=.FALSE.
         END IF
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I8,4X,A,I8)')'ISTATE1=',iRlxRoot,'ISTATE2=',
     &          ISTATE2
           WRITE(6,'(6X,A,I11,4X,A,I11)')'NCI1=',NCI1,'NCI2=',NCI2
         END IF
#ifdef _DEBUGPRINT_
         WRITE(6,*)' TSHinit calls TSHop.'
#endif
         CALL TSHop(CI1,CI2)
#ifdef _DEBUGPRINT_
         WRITE(6,*)' TSHinit back from TSHop.'
#endif
         IF(WFTYP2.EQ.'GENERAL ') THEN
           CALL MkGUGA_Free(SGS(2),CIS(2),EXS(2))
         END IF
         CALL mma_deallocate(CI2)
         CALL mma_deallocate(PART)
      END IF
C
C Check surface hopping to a root higher than the current one
C
      IF (UPROOT) THEN
         ISTATE2=ISTATE1+1
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I3)') 'The upper state is:',ISTATE2
           WRITE(6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',
     &          ENERGY(ISTATE2),' a.u.'
           WRITE(6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',
     &          ENERGY(ISTATE2)-ENERGY(iRlxRoot),' a.u.'
         END IF
*---------------------------------------------------------------------*
* Adapted from RASSI subroutines GTMCTL and READCI                    *
*                                                                     *
C Get wave function parameters for ISTATE2
         JOB2=JBNUM(ISTATE2)
         NACTE2=NACTE(JOB2)
         MPLET2=MLTPLT(JOB2)
         LSYM2=IRREP(JOB2)
         NHOL12=NHOLE1(JOB2)
         NELE32=NELE3(JOB2)
         WFTYP2=RASTYP(JOB2)
         Call NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
         IF(IPGLOB.GE.4) CALL PRPRTTAB(PART)
C For the second wave function
         IF(WFTYP2.EQ.'GENERAL ') THEN
            NRSPRT=3
            DO I=1,8
               NRAS(I,1)=NRS1(I)
               NRAS(I,2)=NRS2(I)
               NRAS(I,3)=NRS3(I)
            END DO
            NRASEL(1)=2*NRS1T-NHOL12
            NRASEL(2)=NACTE2-NELE32
            NRASEL(3)=NACTE2
            CALL SGINIT(NSYM,NACTE2,MPLET2,SGS(2),CIS(2))
            IF(IPGLOB.GT.4) THEN
               WRITE(6,*)'Split-graph structure for JOB2=',JOB2
               CALL SGPRINT(SGS(2))
            END IF
            CALL CXINIT(SGS(2),CIS(2),EXS(2))
C     CI sizes, as function of symmetry, are now known.
            NCI2=CIS(2)%NCSF(LSYM2)
         ELSE
C     Presently, the only other cases are HISPIN, CLOSED or EMPTY.
            NCI2=1
         END IF
         CALL mma_allocate(CI2,NCI2,Label='CI2')
C     Check if the Energy gap is smaller than the threshold.
         Ediff=ABS(ENERGY(ISTATE2)-ENERGY(ISTATE1))
         Ethr=3.0D-2
         IF (Ediff.le.Ethr) THEN
            ChkHop=.TRUE.
         ELSE
            ChkHop=.FALSE.
         END IF
         IF (IPGLOB.GE.2) THEN
           WRITE(6,'(6X,A,I8,4X,A,I8)')'ISTATE1=',iRlxRoot,'ISTATE2=',
     &          ISTATE2
           WRITE(6,'(6X,A,I11,4X,A,I11)')'NCI1=',NCI1,'NCI2=',NCI2
         END IF
#ifdef _DEBUGPRINT_
         WRITE(6,*)' TSHinit calls TSHop.'
#endif
         CALL TSHop(CI1,CI2)
#ifdef _DEBUGPRINT_
         WRITE(6,*)' TSHinit back from TSHop.'
#endif
         IF(WFTYP2.EQ.'GENERAL ') THEN
           CALL MkGUGA_Free(SGS(2),CIS(2),EXS(2))
         END IF
         CALL mma_deallocate(CI2)
         CALL mma_deallocate(PART)
      END IF
*
      IF(WFTYP1.EQ.'GENERAL ') THEN
        CALL MkGUGA_Free(SGS(1),CIS(1),EXS(1))
      END IF
      CALL mma_deallocate(CI1)
*
      END SUBROUTINE TSHinit
