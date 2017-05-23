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
      SUBROUTINE TSHinit()
C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='TSHINIT')
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "tshcntrl.fh"
      DIMENSION ISGSTR1(NSGSIZE), ISGSTR2(NSGSIZE)
      DIMENSION ICISTR1(NCISIZE), ICISTR2(NCISIZE)
      DIMENSION IXSTR1(NXSIZE), IXSTR2(NXSIZE)
      INTEGER      I,JOB1,JOB2,iRlxRoot
      CHARACTER*8  WFTYP1,WFTYP2
      LOGICAL   LOWROOT, UPROOT
*


      CALL QENTER(ROUTINE)
C
C Print a banner
C
      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(6X,100A1)') ('*',i=1,100)
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A,36X,A,37X,A)')
     &       '*',' Surface hopping section ','*'
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,100A1)') ('*',i=1,100)
        WRITE(6,*)
        WRITE(6,*)

        WRITE(6,'(6X,A)')'Surface hopping section'
        WRITE(6,'(6X,A)')'-----------------------'
      END IF
C
C Get the current state and print it's energy
C
      CALL Get_iScalar('Relax CASSCF root',iRlxRoot)
      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,'(6X,A,I14)')'The current state is:',iRlxRoot
        WRITE(6,'(6X,A,6X,E15.6,A,/)')'Its energy is:',ENERGY(iRlxRoot),
     &                              ' a.u.'
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
          CALL SGINIT(NSYM,NACTE1,MPLET1,NRSPRT,NRAS,NRASEL,ISGSTR1)
          IF(IPGLOB.GT.DEBUG) THEN
             WRITE(6,*)'Split-graph structure for JOB1=',JOB1
             CALL SGPRINT(ISGSTR1)
          END IF
          CALL SGSVAL(ISGSTR1,NSYM,NASHT,LISM,NVERT,LDRT,
     &               LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
          CALL CXINIT(ISGSTR1,ICISTR1,IXSTR1)
          CALL CXSVAL(ICISTR1,IXSTR1,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &               LNOCSF,LIOCSF,NWALK,LICASE,
     &               MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &               LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C CI sizes, as function of symmetry, are now known.
          NCI1=IWORK(LNCSF-1+LSYM1)
      ELSE
* The only other cases are HISPIN, CLOSED or EMPTY.
* NOTE: The HISPIN case is suspected to be buggy. Not used now.
          NCI1=1
      END IF
      CALL GETMEM('GTDMCI1','ALLO','REAL',LCI1,NCI1)

C
C Check for the possibility to hop to a lower root
C
      IF ((ISTATE1-1).GE.1) THEN
         LOWROOT=.TRUE.
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I2)') "There is a lower root, which is: "
     &     ,LROOT(ISTATE1-1)
         END IF
      ELSE
         LOWROOT=.FALSE.
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A)') "There is no lower root"
         END IF
      END IF
C
C Check for the possibility to hop to an upper root
C
      IF ((ISTATE1+1).LE.NSTATE) THEN
         UPROOT=.TRUE.
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I2,/)') "There is an upper root, which is: "
     &     ,LROOT(ISTATE1+1)
         END IF
      ELSE
         UPROOT=.FALSE.
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,/)') "There is no upper root"
         END IF
      END IF
C
C Check surface hopping to a root lower than the current one
C
      IF (LOWROOT) THEN
         ISTATE2=ISTATE1-1
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I3)') 'The lower state is:',ISTATE2
           WRITE(6,'(6X,A,6X,E15.6,A)') 'Its energy is:',
     &          ENERGY(ISTATE2),' a.u.'
           WRITE(6,'(6X,A,12X,E15.6,A,/)') 'Ediff = ',
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
         LPART=NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
         IF(IPGLOB.GE.DEBUG) CALL PRPRTTAB(IWORK(LPART))
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
            CALL SGINIT(NSYM,NACTE2,MPLET2,NRSPRT,NRAS,
     &           NRASEL,ISGSTR2)
            IF(IPGLOB.GT.DEBUG) THEN
               WRITE(6,*)'Split-graph structure for JOB2=',JOB2
               CALL SGPRINT(ISGSTR2)
            END IF
            CALL SGSVAL(ISGSTR2,NSYM,NASHT,LISM,NVERT,LDRT,
     &           LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
            CALL CXINIT(ISGSTR2,ICISTR2,IXSTR2)
            CALL CXSVAL(ICISTR2,IXSTR2,NMIDV,NIPWLK,LNOW,LIOW,
     &           LNCSF,LNOCSF,LIOCSF,NWALK,LICASE,MXEO,
     &           LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,LVTAB,
     &           LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C     CI sizes, as function of symmetry, are now known.
            NCI2=IWORK(LNCSF-1+LSYM2)
         ELSE
C     Presently, the only other cases are HISPIN, CLOSED or EMPTY.
            NCI2=1
         END IF
         CALL GETMEM('GTDMCI2','ALLO','REAL',LCI2,NCI2)
C     Check if the Energy gap is smaller than the threshold.
         Ediff=ABS(ENERGY(ISTATE2)-ENERGY(ISTATE1))
         Ethr=3.0D-2
         IF (Ediff.le.Ethr) THEN
            ChkHop=.TRUE.
         ELSE
            ChkHop=.FALSE.
         END IF
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I8,4X,A,I8)')'ISTATE1=',iRlxRoot,'ISTATE2=',
     &          ISTATE2
           WRITE(6,'(6X,A,I11,4X,A,I11)')'NCI1=',NCI1,'NCI2=',NCI2
         END IF
#ifdef _DEBUG_
         WRITE(6,*)' TSHinit calls TSHop.'
#endif
         CALL TSHop(WORK(LCI1),WORK(LCI2))
#ifdef _DEBUG_
         WRITE(6,*)' TSHinit back from TSHop.'
#endif
         IF(WFTYP2.EQ.'GENERAL ') THEN
           CALL CXCLOSE(ISGSTR2,ICISTR2,IXSTR2)
           CALL SGCLOSE(ISGSTR2)
         END IF
         CALL GETMEM('GTDMCI2','FREE','REAL',LCI2,NCI2)
         CALL KILLOBJ(LPART)
      END IF
C
C Check surface hopping to a root higher than the current one
C
      IF (UPROOT) THEN
         ISTATE2=ISTATE1+1
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I3)') 'The upper state is:',ISTATE2
           WRITE(6,'(6X,A,6X,E15.6,A)') 'Its energy is:',
     &          ENERGY(ISTATE2),' a.u.'
           WRITE(6,'(6X,A,12X,E15.6,A,/)') 'Ediff = ',
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
         LPART=NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
         IF(IPGLOB.GE.DEBUG) CALL PRPRTTAB(IWORK(LPART))
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
            CALL SGINIT(NSYM,NACTE2,MPLET2,NRSPRT,NRAS,
     &           NRASEL,ISGSTR2)
            IF(IPGLOB.GT.DEBUG) THEN
               WRITE(6,*)'Split-graph structure for JOB2=',JOB2
               CALL SGPRINT(ISGSTR2)
            END IF
            CALL SGSVAL(ISGSTR2,NSYM,NASHT,LISM,NVERT,LDRT,
     &           LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
            CALL CXINIT(ISGSTR2,ICISTR2,IXSTR2)
            CALL CXSVAL(ICISTR2,IXSTR2,NMIDV,NIPWLK,LNOW,LIOW,
     &           LNCSF,LNOCSF,LIOCSF,NWALK,LICASE,MXEO,
     &           LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,LVTAB,
     &           LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C     CI sizes, as function of symmetry, are now known.
            NCI2=IWORK(LNCSF-1+LSYM2)
         ELSE
C     Presently, the only other cases are HISPIN, CLOSED or EMPTY.
            NCI2=1
         END IF
         CALL GETMEM('GTDMCI2','ALLO','REAL',LCI2,NCI2)
C     Check if the Energy gap is smaller than the threshold.
         Ediff=ABS(ENERGY(ISTATE2)-ENERGY(ISTATE1))
         Ethr=3.0D-2
         IF (Ediff.le.Ethr) THEN
            ChkHop=.TRUE.
         ELSE
            ChkHop=.FALSE.
         END IF
         IF (IPGLOB.GE.USUAL) THEN
           WRITE(6,'(6X,A,I8,4X,A,I8)')'ISTATE1=',iRlxRoot,'ISTATE2=',
     &          ISTATE2
           WRITE(6,'(6X,A,I11,4X,A,I11)')'NCI1=',NCI1,'NCI2=',NCI2
         END IF
#ifdef _DEBUG_
         WRITE(6,*)' TSHinit calls TSHop.'
#endif
         CALL TSHop(WORK(LCI1),WORK(LCI2))
#ifdef _DEBUG_
         WRITE(6,*)' TSHinit back from TSHop.'
#endif
         IF(WFTYP2.EQ.'GENERAL ') THEN
           CALL CXCLOSE(ISGSTR2,ICISTR2,IXSTR2)
           CALL SGCLOSE(ISGSTR2)
         END IF
         CALL GETMEM('GTDMCI2','FREE','REAL',LCI2,NCI2)
         CALL KILLOBJ(LPART)
      END IF
*
C      ELSE
C         LCI2=LCI1
C         NCI2=NCI1
C         CALL GETMEM('GTDMCI2','ALLO','REAL',LCI2,NCI2)
C      END IF
      IF(WFTYP1.EQ.'GENERAL ') THEN
        CALL CXCLOSE(ISGSTR1,ICISTR1,IXSTR1)
        CALL SGCLOSE(ISGSTR1)
      END IF
      CALL GETMEM('GTDMCI1','FREE','REAL',LCI1,NCI1)
      CALL QEXIT(ROUTINE)
      RETURN
*
      END
