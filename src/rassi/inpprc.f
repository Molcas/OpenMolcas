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
      SUBROUTINE INPPRC
      use rassi_global_arrays, only: HAM, ESHFT
      use rassi_aux, Only : jDisk_TDM, AO_Mode, JOB_INDEX, CMO1, CMO2,
     &                      DMAB, mTRA
      use kVectors
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='INPPRC')
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "rasdim.fh"
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "SysDef.fh"
#include "centra.fh"
      CHARACTER*8 LABEL
      CHARACTER*8 LABEL2
      !CHARACTER*8 LABEL3
      CHARACTER*8 PRPLST(MXPROP)
      Character*3 lIrrep(8)
      INTEGER ICMPLST(MXPROP)
      LOGICAL JOBMATCH
      DIMENSION IDUM(1)
* Analysing and post-processing the input that was read in readin_rassi.

      CALL QENTER(ROUTINE)

      Call mma_allocate(jDisk_TDM,2,nState*(nState+1)/2,
     &                  Label='jDisk_TDM')
      jDisk_TDM(:,:)=-1
* PAM07: The printing of spin-orbit Hamiltonian matrix elements:
* If no value for SOTHR_PRT was given in the input, it has a
* negative value that was set in init_rassi:
      IF(SOTHR_PRT.lt.0.0D0) THEN
* Assign default settings. SOTHR_PRT is in cm-1:
       IF(IPGLOB.ge.DEBUG) THEN
        NSOTHR_PRT=10000
        SOTHR_PRT=0.0001D0
      ELSE IF (IPGLOB.ge.VERBOSE) THEN
        NSOTHR_PRT=100
        SOTHR_PRT=0.01D0
      ELSE IF(IPGLOB.ge.USUAL) THEN
        NSOTHR_PRT=20
        SOTHR_PRT=1.00D0
       ELSE
        NSOTHR_PRT=0
       END IF
      END IF

C Some sizes:
      NBSQ=0
      NBMX=0
      IPRP=0
      DO I=1,NSYM
        NBI=NBASF(I)
        NBMX=MAX(NBMX,NBI)
        NBSQPR(I)=NBSQ
        NBSQ=NBSQ+NBI**2
      END DO
      NBTRI=(NBSQ+NBST)/2
      NLEV=0
      DO I=1,NSYM
        NLEV=NLEV+NRS1(I)+NRS2(I)+NRS3(I)
      END DO
C Sizes of some data sets:
C ACTUAL SIZES OF TDMAB AND TDMZZ DEPENDS ON BOTH LSYM1 AND LSYM2.
C HOWEVER, MAX POSSIBLE SIZE IS WHEN LSYM1=LSYM2.
      NCMO=NOSH(1)*NBASF(1)
      NTRA=NOSH(1)**2
      NTDMZZ=NBASF(1)**2
      DO IS=2,NSYM
        NCMO=NCMO+NOSH(IS)*NBASF(IS)
        NTRA=NTRA+NOSH(IS)**2
        NTDMZZ=NTDMZZ+NBASF(IS)**2
      END DO
      NTDMS=(NTDMZZ+NBST)/2
      NTDMA=NTDMS
      NTDMAB=NTRA
      IF((SONATNSTATE.GT.0).OR.NATO.OR.Do_TMOM) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info: creating TDMFILE'
        LUTDM=21
        LUTDM=IsFreeUnit(LUTDM)
        FNTDM='TDMFILE'
        CALL DANAME_MF(LUTDM,FNTDM)
        AO_Mode=.True.
        iByte=8*3*nstate*(nstate-1)/2*nTDMZZ
*
*       For the time we will move over to compact mode if the required
*       estimate of disk space if more than 1 Gb.
*
        If (iByte.gt.1024**3) AO_Mode=.False.
*       Force for debugging purpose.
        If (Force_NON_AO_TDM)  AO_Mode=.False.
        WRITE(6,*) '       estimated file size ', iByte/1024, 'kB'
*
*       For small basis set with symmetry we might not benefit from
*       doing this.
*
        If (NASHT**2+1.gt.nTDMAB) AO_Mode=.True.
*
        If (.NOT.AO_Mode) Then
           WRITE(6,*) '       TDMs in reduced format'
           Call mma_allocate(JOB_INDEX,nState,Label='JOB_INDEX')
           Call ICopy(nState,iWork(lJBNUM),1,JOB_INDEX,1)
*          Write (6,*) 'Job_Index=',Job_Index
           Call mma_allocate(CMO1,nCMO,Label='CMO1')
           Call mma_allocate(CMO2,nCMO,Label='CMO2')
           Call mma_allocate(DMAB,nTDMAB,Label='DMAB')
           mTRA=nTRA
        Else
           WRITE(6,*) '       TDMs in AO format'
        End If
        WRITE(6,*)
      END IF

C Upcase property names in lists of requests:
      DO IPROP=1,NPROP
        CALL UPCASE(PNAME(IPROP))
      END DO
      DO ISOPR=1,NSOPR
        CALL UPCASE(SOPRNM(ISOPR))
      END DO

C Which properties are available in the ONEINT file?
C (IPUSED will be set later, set it to zero now.)
      !write(6,*)"Which properties are available in the ONEINT file?"
      IPRP=0
      IRC=-1
      IOPT=15
      LABEL='UNDEF'
      CALL iRDONE(IRC,IOPT,LABEL,ICMP,IDUM,ISYLAB)
      IF(IRC.EQ.0) NSIZ=IDUM(1)
      IF(IRC.NE.0) GOTO 110
      IPRP=1
      CALL UPCASE(LABEL)
      PRPLST(1)=LABEL
      ICMPLST(1)=ICMP
      IPUSED(1)=0
*
      DO I=1,MXPROP
        IF(IPRP.GE.MXPROP) GOTO 110
        IRC=-1
        IOPT=23
        CALL iRDONE(IRC,IOPT,LABEL,ICMP,IDUM,ISYLAB)
        IF(IRC.EQ.0) NSIZ=IDUM(1)
        IF(IRC.NE.0) GOTO 110
        IPRP=IPRP+1
        CALL UPCASE(LABEL)
        PRPLST(IPRP)=LABEL
        ICMPLST(IPRP)=ICMP
        IPUSED(IPRP)=0

c Copy the EF2 integral label for hyperfine calculations
        IF(LABEL(1:3).EQ.'EF2') THEN
          IF(IPRP.GE.MXPROP) GOTO 110
          IPRP=IPRP+1
          LABEL2=LABEL
          LABEL2(1:3)='ASD'
          PRPLST(IPRP)=LABEL2
          ICMPLST(IPRP)=ICMP
          IPUSED(IPRP)=0

        END IF

        IF(LABEL(1:4).EQ.'PSOI') THEN
          IF(IPRP.GE.MXPROP) GOTO 110
          IPRP=IPRP+1
          LABEL2=LABEL
          LABEL2(1:4)='PSOP'
          PRPLST(IPRP)=LABEL2
          ICMPLST(IPRP)=ICMP
          IPUSED(IPRP)=0
        END IF
        IF(LABEL(1:6).EQ.'DMS  1') THEN
          IF(IPRP.GE.MXPROP) GOTO 110
          IPRP=IPRP+1
          LABEL2=LABEL
          LABEL2(1:6)='DMP   '
          PRPLST(IPRP)=LABEL2
          ICMPLST(IPRP)=ICMP
          IPUSED(IPRP)=0
        END IF
      END DO
110   CONTINUE
      NPRPLST=IPRP
*
*     Add empty slots for on-the-fly TM integrals.
*
*     If the RASSI code is run several instances on the same job some
*     of these labels will already be available on the file and need
*     not to be added to the list.
*
      IF (Do_TMOM.AND.PRPLST(IPRP)(1:4).NE.'TMOM') THEN
         PRPLST(IPRP+ 1)='TMOM0  R'
         ICMPLST(IPRP+ 1)=1
         IPUSED(IPRP+ 1)=0
         PRPLST(IPRP+ 2)='TMOM0  I'
         ICMPLST(IPRP+ 2)=1
         IPUSED(IPRP+ 2)=0
         IPRP=IPRP+2
*
         PRPLST(IPRP+ 1)='TMOM  RS'
         ICMPLST(IPRP+ 1)=1
         IPUSED(IPRP+ 1)=0
         PRPLST(IPRP+ 2)='TMOM  RS'
         ICMPLST(IPRP+ 2)=2
         IPUSED(IPRP+ 2)=0
         PRPLST(IPRP+ 3)='TMOM  RS'
         ICMPLST(IPRP+ 3)=3
         IPUSED(IPRP+ 3)=0
         PRPLST(IPRP+ 4)='TMOM  RA'
         ICMPLST(IPRP+ 4)=1
         IPUSED(IPRP+ 4)=0
         PRPLST(IPRP+ 5)='TMOM  RA'
         ICMPLST(IPRP+ 5)=2
         IPUSED(IPRP+ 5)=0
         PRPLST(IPRP+ 6)='TMOM  RA'
         ICMPLST(IPRP+ 6)=3
         IPUSED(IPRP+ 6)=0
         PRPLST(IPRP+ 7)='TMOM  IS'
         ICMPLST(IPRP+ 7)=1
         IPUSED(IPRP+ 7)=0
         PRPLST(IPRP+ 8)='TMOM  IS'
         ICMPLST(IPRP+ 8)=2
         IPUSED(IPRP+ 8)=0
         PRPLST(IPRP+ 9)='TMOM  IS'
         ICMPLST(IPRP+ 9)=3
         IPUSED(IPRP+ 9)=0
         PRPLST(IPRP+10)='TMOM  IA'
         ICMPLST(IPRP+10)=1
         IPUSED(IPRP+10)=0
         PRPLST(IPRP+11)='TMOM  IA'
         ICMPLST(IPRP+11)=2
         IPUSED(IPRP+11)=0
         PRPLST(IPRP+12)='TMOM  IA'
         ICMPLST(IPRP+12)=3
         IPUSED(IPRP+12)=0
         IPRP=IPRP+12
*
C Not currently in use.
C        PRPLST(IPRP+ 1)='TMOM2  R'
C        ICMPLST(IPRP+ 1)=1
C        IPUSED(IPRP+ 1)=0
C        CALL MKTDAB(0.0D0,WERD,WDMAB,iRC)PRPLST(IPRP+ 2)='TMOM2  R'
C        ICMPLST(IPRP+ 2)=2
C        IPUSED(IPRP+ 2)=0
C        PRPLST(IPRP+ 3)='TMOM2  R'
C        ICMPLST(IPRP+ 3)=3
C        IPUSED(IPRP+ 3)=0
C        PRPLST(IPRP+ 4)='TMOM2  I'
C        ICMPLST(IPRP+ 4)=1
C        IPUSED(IPRP+ 4)=0
C        PRPLST(IPRP+ 5)='TMOM2  I'
C        ICMPLST(IPRP+ 5)=2
C        IPUSED(IPRP+ 5)=0
C        PRPLST(IPRP+ 6)='TMOM2  I'
C        ICMPLST(IPRP+ 6)=3
C        IPUSED(IPRP+ 6)=0
C        IPRP=IPRP+6
C     Else IF (Do_TMOM) Then
C        PRPLST(IPRP+ 1)='TMOM0  R'
C        ICMPLST(IPRP+ 1)=1
C        IPUSED(IPRP+ 1)=0
C        PRPLST(IPRP+ 2)='TMOM0  I'
C        ICMPLST(IPRP+ 2)=1
C        IPUSED(IPRP+ 2)=0
C        IPRP=IPRP+2
*
C Not currently in use.
C        PRPLST(IPRP+ 1)='TMOM2  R'
C        ICMPLST(IPRP+ 1)=2
C        IPUSED(IPRP+ 1)=0
C        PRPLST(IPRP+ 2)='TMOM2  R'
C        ICMPLST(IPRP+ 2)=3
C        IPUSED(IPRP+ 2)=0
C        PRPLST(IPRP+ 3)='TMOM2  I'
C        ICMPLST(IPRP+ 3)=2
C        IPUSED(IPRP+ 3)=0
C        PRPLST(IPRP+ 4)='TMOM2  I'
C        ICMPLST(IPRP+ 4)=3
C        IPUSED(IPRP+ 4)=0
C        IPRP=IPRP+4
      END IF
*
      NPRPLST=IPRP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C Add some property names by defaults, if no input:                   C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF (NPROP.EQ.0) THEN
C
         IF (NSOPR.EQ.0) THEN
C If no input at all, use this selection:
            DO IPRP=1,NPRPLST
C
               IF (PRPLST(IPRP).eq.'MLTPL  0' .or.
     &             PRPLST(IPRP).eq.'MLTPL  1' .or.
     &             PRPLST(IPRP).eq.'MLTPL  2' .or.
     &             PRPLST(IPRP).eq.'MLTPL  3' .or.
     &             PRPLST(IPRP).eq.'OMQ     ' .or.
     &             PRPLST(IPRP).eq.'ANGMOM' .or.
     &             PRPLST(IPRP)(1:4).eq.'TMOM' .or.
     &             PRPLST(IPRP).eq.'VELOCITY' .or.
     &             PRPLST(IPRP).eq.'MLTPV  2' .or.
     &             PRPLST(IPRP)(1:4).eq.'EMFR') THEN
C
                  NSOPR=NSOPR+1
                  SOPRNM(NSOPR)=PRPLST(IPRP)
C
C SET the proper default type label. That is, which WE-reduced density
C should the integral be contracted with. If not set here it will
C default to a Hermitian singlet reduced density.
                  IF (PRPLST(IPRP).EQ.'ANGMOM'  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'VELOCITY'  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'OMQ     '  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'MLTPV  2'  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM  RS'  ) THEN
                     SOPRTP(NSOPR)='HERMSING'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM  RA'  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM  IS'  ) THEN
                     SOPRTP(NSOPR)='HERMSING'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM  IA'  ) THEN
                     SOPRTP(NSOPR)='ANTISING'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM0  R'  ) THEN
                     SOPRTP(NSOPR)='HERMTRIP'
                  ELSE IF (PRPLST(IPRP).EQ.'TMOM0  I'  ) THEN
                     SOPRTP(NSOPR)='ANTITRIP'
                  ELSE
                     SOPRTP(NSOPR)='HERMSING'
                  END IF
C
C Set the number of elements of the property
                  ISOCMP(NSOPR)=ICMPLST(IPRP)
C
C Add properties for the explicit spin part of the transition moments
C in the case of spin-orbit coupled wave functrions.
C
                  IF (IFSO) THEN
                     IF (PRPLST(IPRP).EQ.'MLTPL  0') THEN
                        NSOPR=NSOPR+1
                        SOPRNM(NSOPR)=PRPLST(IPRP)
                        ISOCMP(NSOPR)=ICMPLST(IPRP)
                        SOPRTP(NSOPR)='ANTITRIP'
                     ELSE IF (PRPLST(IPRP).EQ.'MLTPL  1') THEN
                        NSOPR=NSOPR+1
                        SOPRNM(NSOPR)=PRPLST(IPRP)
                        ISOCMP(NSOPR)=ICMPLST(IPRP)
                        SOPRTP(NSOPR)='ANTITRIP'
C Uncomment to activate more options!!!
C                    ELSE IF (PRPLST(IPRP).EQ.'MLTPL  2') THEN
C                       NSOPR=NSOPR+1
C                       SOPRNM(NSOPR)=PRPLST(IPRP)
C                       ISOCMP(NSOPR)=ICMPLST(IPRP)
C                       SOPRTP(NSOPR)='HERMTRIP'
C                    ELSE IF (PRPLST(IPRP).EQ.'ANGMOM'  ) THEN
C                       NSOPR=NSOPR+1
C                       SOPRNM(NSOPR)=PRPLST(IPRP)
C                       ISOCMP(NSOPR)=ICMPLST(IPRP)
C                       SOPRTP(NSOPR)='ANTITRIP'
                     END IF
                  END IF
C
               END IF
C
C Add some properties if DQVD is requested
               IF (DQVD) THEN
C                 'MLTPL  2' already there by default
                  IF (PRPLST(IPRP).eq.'EF0    1') THEN
                     NSOPR=NSOPR+1
                     SOPRNM(NSOPR)=PRPLST(IPRP)
                     ISOCMP(NSOPR)=ICMPLST(IPRP)
                  END IF
               END IF
C
            END DO
         END IF
C If no PROP input, copy the SOPR selection:
         NPROP=NSOPR
         DO IPROP=1,NPROP
            PNAME(IPROP)=SOPRNM(IPROP)
            PTYPE(IPROP)=SOPRTP(IPROP)
            ICOMP(IPROP)=ISOCMP(IPROP)
         END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
      ELSE
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C If no SOPR input, copy the PROP selection:
         IF (NSOPR.EQ.0) THEN
            NSOPR=NPROP
            DO ISOPR=1,NSOPR
               SOPRNM(ISOPR)=PNAME(ISOPR)
               SOPRTP(ISOPR)=PTYPE(ISOPR)
               ISOCMP(ISOPR)=ICOMP(ISOPR)
            END DO
         END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
      END IF
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
* Lists (above) now contain either a default choice, or
* a selection by the user. Check that integrals are
* available on the oneint file.

* Check if we need to activate IFJ2/IFJZ automatically
      if(natoms.eq.1) then
       ifj2=1
      endif
      xaxis=0.d0
      zaxis=0.d0
      Do iatom=1,natoms
       xaxis=xaxis+abs(coor(1,iatom))+abs(coor(2,iatom))
       zaxis=zaxis+abs(coor(3,iatom))
      Enddo
      if(xaxis.lt.1.d-10.and.zaxis.gt.0.1d0) ifjz=1

* Check if angular momentum integrals have been computed
      MISSAMX=1
      MISSAMY=1
      MISSAMZ=1
      DO IPRP=1,NPRPLST
         IF (PRPLST(IPRP).eq.'ANGMOM  ') THEN
            IF (ICMPLST(IPRP).eq.1) MISSAMX=0
            IF (ICMPLST(IPRP).eq.2) MISSAMY=0
            IF (ICMPLST(IPRP).eq.3) MISSAMZ=0
         END IF
      END DO
      IF (MISSAMX+MISSAMY+MISSAMZ.gt.0) THEN
       IF(IFJ2.eq.1) THEN
        write(6,*)' J2 values cannot be computed.'
        write(6,*)' Reason: Angular momentum integrals are missing.'
        IFJ2=0
       END IF
       IF(IFGCAL .OR. IFXCAL) THEN
        write(6,*)' Neither GCAL or XCAL can be computed.'
        write(6,*)' Reason: Angular momentum integrals are missing.'
        IFJ2=0
       END IF
      END IF
      IF (MISSAMZ.gt.0) THEN
       IF(IFJZ.eq.1) THEN
        write(6,*)' Omega values will not be computed.'
        write(6,*)' Reason: Angular momentum integrals are missing.'
        IFJZ=0
       END IF
      END IF

C Is everything available that we may need?
      IMiss=0
      DO IPROP=1,NPROP
       DO IPRP=1,NPRPLST
        IF(PNAME(IPROP).EQ.PRPLST(IPRP).AND.
     &        ICOMP(IPROP).EQ.ICMPLST(IPRP)) THEN
         IPUSED(IPRP)=1
         GOTO 120
        END IF
       END DO

       IMiss=IMiss+1
       IF(IMiss.eq.1) THEN
         Write(6,*)
         Call WarningMessage(1,'Requested integrals are missing.')
         Write(6,*)' Property name, and component:',
     &           PNAME(IPROP),ICOMP(IPROP)
         Write(6,*)' This record cannot be found. Some of the requested'
         Write(6,*)' properties cannot be computed. Suggested fix: Try'
         Write(6,*)' recomputing one-electron integrals with keyword'
         Write(6,*)' ''OneOnly'', and additional keywords for the'
         Write(6,*)' properties needed.'
       ELSE
         Write(6,*)' Also missing:',PNAME(IPROP),ICOMP(IPROP)
       END IF

 120   CONTINUE
      END DO
      IMiss=0
      DO ISOPR=1,NSOPR
       DO IPRP=1,NPRPLST
        IF(SOPRNM(ISOPR).EQ.PRPLST(IPRP).AND.
     &        ISOCMP(ISOPR).EQ.ICMPLST(IPRP)) THEN
         IPUSED(IPRP)=1
         GOTO 130
        END IF
       END DO

       IMiss=IMiss+1
       IF(IMiss.eq.1) THEN
         Write(6,*)
         Call WarningMessage(1,'Requested integrals are missing.')
         Write(6,*)' SO-Property name, and component:',
     &           SOPRNM(ISOPR),ISOCMP(ISOPR)
         Write(6,*)' This record cannot be found. Some of the requested'
         Write(6,*)' properties cannot be computed. Suggested fix: Try'
         Write(6,*)' recomputing one-electron integrals with keyword'
         Write(6,*)' ''OneOnly'', and additional keywords for the'
         Write(6,*)' properties needed.'
       ELSE
         Write(6,*)' Also missing:',SOPRNM(ISOPR),ISOCMP(ISOPR)
       END IF

 130   CONTINUE
      END DO
cnf
      If (IfDCpl) Then
         Call Get_iScalar('Unique atoms',natom)
         IErr = 3*natom
         Do IPrp = 1, NPrpLst
            If(PrpLst(IPrp)(1:3).eq.'EF1') Then
               IErr=IErr -1
               IPUsed(IPrp)=1
            End If
         End Do
         If (IErr.ne.0) Then
            Write(6,*)
            Write(6,'(A,i5,A)') '  Approx derivative couplings require',
     &                 3*natom,' field integrals.'
            Write(6,*)' Add the keywords EFLD=0 and ONEONLY to'//
     &                ' the SEWARD input and recompute ONEINT.'
            Write(6,*)' Program will continue, but DerCpl calculations'
     &                //' are disabled.'
            IfDCpl = .False.
         End If
      End If
cnf
C Temporary use IPUSED to mark operators that may be needed:
      DO IPRP=1,NPRPLST
       DO IPROP=1,NPROP
        IF(PNAME(IPROP).EQ.PRPLST(IPRP).AND.
     &        ICOMP(IPROP).EQ.ICMPLST(IPRP)) THEN
          IPUSED(IPRP)=1
          GOTO 221
        END IF
       END DO
 221   CONTINUE
      END DO

      IF(.not.IFHAM) THEN
        IF(IFSHFT) THEN
          Call WarningMessage(1,'Ignored user request.')
          WRITE(6,*)' INPCTL Warning: Both keywords ONEL and SHIFT.'
          WRITE(6,*)' User-supplied diagonal energy shifts are'//
     &              ' meaningless since no Hamiltonian will be'
          WRITE(6,*)' used/computed anyway. Ignored.'
          IFSHFT=.FALSE.
        END IF
        IF(IFHDIA) THEN
          Call WarningMessage(1,'Ignored user request.')
          WRITE(6,*)' INPCTL Warning: Both keywords ONEL and HDIAG.'
          WRITE(6,*)' User-supplied diagonal H elements are'//
     &              ' meaningless since no Hamiltonian will be'
          WRITE(6,*)' used/computed anyway. Ignored.'
          IFHDIA=.FALSE.
        END IF
      END IF
C In any Spin-Orbit calculation, we need SO integrals as well.
C Right now, only AMFI is available. Check that first.
      IF(IFSO.OR.IFGCAL.OR.IFXCAL.OR.IFMCAL) THEN
        IERR=3
        DO IPRP=1,NPRPLST
          IF(PRPLST(IPRP).EQ.'AMFI    ') THEN
           IERR=IERR-1
           IPUSED(IPRP)=1

          END IF
          !IF(PRPLST(IPRP).EQ.'PSOI    ') THEN
          ! IERR=IERR-1
          ! IPUSED(IPRP)=1
          !END IF

        END DO
        IF(IERR.NE.0) THEN
          Call WarningMessage(1,'Incomplete integrals.')
          WRITE(6,*)' Spin-Orbit interaction calculation was requested'
          WRITE(6,*)' but this requires three components of AMFI'
          WRITE(6,*)' integrals. Add the keywords AMFI and ONEONLY to'
          WRITE(6,*)' the SEWARD input and recompute ONEINT.'
          WRITE(6,*)' Program will continue, but SO/EPRG/MAGN'
          WRITE(6,*)' calculations are disabled.'
          IFSO=.FALSE.
          NSOPR=0
          IFGCAL=.FALSE.
          IFXCAL=.FALSE.
          IFMCAL=.FALSE.
        END IF
      END IF
C SVC2009 Check for the presence of the ANGMOM integrals needed for G factor
C or Magnetic Moment calculations.
      IF(IFGCAL.OR.IFXCAL.OR.IFMCAL) THEN
        IERR=3
        DO IPRP=1,NPRPLST
          IF(PRPLST(IPRP).EQ.'ANGMOM  ') THEN
           IERR=IERR-1
           IPUSED(IPRP)=1
          END IF

          !IF(PRPLST(IPRP).EQ.'PSOI    ') THEN
          ! write(6,*)"5*****rassi/inpprc ANGMOM"
          ! IERR=IERR-1
          ! IPUSED(IPRP)=1
          !END IF

        END DO
        IF(IERR.NE.0) THEN
          Call WarningMessage(1,'Incomplete integrals.')
          WRITE(6,*)' EPRG or MAGN keyword was requested'
          WRITE(6,*)' but this requires three components of ANGMOM'
          WRITE(6,*)' integrals. Add the keywords ANGMOM and ONEONLY to'
          WRITE(6,*)' the SEWARD input and recompute ONEINT.'
          WRITE(6,*)' Program will continue, but EPRG/MAGN calculations'
          WRITE(6,*)' are disabled.'
          IFGCAL=.FALSE.
          IFXCAL=.FALSE.
          IFMCAL=.FALSE.
        END IF
      END IF
C Similarly, check integrals for which we want matrix elements over
C SO eigenstates.
      DO IPRP=1,NPRPLST
       DO ISOPR=1,NSOPR
        IF(SOPRNM(ISOPR).EQ.PRPLST(IPRP).AND.
     &        ISOCMP(ISOPR).EQ.ICMPLST(IPRP)) THEN
          IPUSED(IPRP)=1
          GOTO 222
        END IF
       END DO
 222   CONTINUE
      END DO
C
C Reassemble the PNAME, ICOMP arrays.
C
      DO IPRP=1,NPRPLST
         IF (IPUSED(IPRP).EQ.0) THEN
            DO IPROP = 1, NPROP
               IF (PNAME(IPROP).EQ.PRPLST(IPRP)) PNAME(IPROP)='REMOVE'
            END DO
         ELSE
            IADD=1
            DO IPROP = 1, NPROP
               IF (PNAME(IPROP).EQ.PRPLST(IPRP) .AND.
     &             ICOMP(IPROP).EQ.ICMPLST(IPRP)) IADD=0
            END DO
            IF (IADD.EQ.1) THEN
               NPROP=NPROP+1
               PNAME(NPROP)=PRPLST(IPRP)
               PTYPE(NPROP)='UNDEF.  '
               ICOMP(NPROP)=ICMPLST(IPRP)
            END IF
         END IF
      END DO
      MPROP=NPROP
      DO IPROP = 1, NPROP
         IF (PNAME(IPROP).EQ.'DONE') EXIT
         IF (PNAME(IPROP).EQ.'REMOVE') THEN
            DO JPROP = IPROP, MPROP-1
               PNAME(JPROP)=PNAME(JPROP+1)
               PTYPE(JPROP)=PTYPE(JPROP+1)
               ICOMP(JPROP)=ICOMP(JPROP+1)
            END DO
            PNAME(MPROP)='DONE'
            MPROP=MPROP-1
         END IF
      END DO
      NPROP=MPROP
C
C Reassemble the SOPRNM, ISOCMP arrays:
C
      DO IPRP=1,NPRPLST
         IF (IPUSED(IPRP).EQ.0) THEN
            DO ISOPR = 1, NSOPR
               IF (SOPRNM(ISOPR).EQ.PRPLST(IPRP)) SOPRNM(ISOPR)='REMOVE'
            END DO
         ELSE
            IADD=1
            DO ISOPR = 1, NSOPR
               IF (SOPRNM(ISOPR).EQ.PRPLST(IPRP) .AND.
     &             ISOCMP(ISOPR).EQ.ICMPLST(IPRP)) IADD=0
            END DO
            IF (IADD.EQ.1) THEN
               NSOPR=NSOPR+1
               SOPRNM(NSOPR)=PRPLST(IPRP)
               SOPRTP(NSOPR)='UNDEF.  '
               ISOCMP(NSOPR)=ICMPLST(IPRP)
            END IF
         END IF
      END DO
      MSOPR=NSOPR
      DO ISOPR = 1, NSOPR
         IF (SOPRNM(ISOPR).EQ.'DONE') EXIT
         IF (SOPRNM(ISOPR).EQ.'REMOVE') THEN
            DO JSOPR = ISOPR, MSOPR-1
               SOPRNM(JSOPR)=SOPRNM(JSOPR+1)
               SOPRTP(JSOPR)=SOPRTP(JSOPR+1)
               ISOCMP(JSOPR)=ISOCMP(JSOPR+1)
            END DO
            PNAME(MSOPR)='DONE'
            MSOPR=MSOPR-1
         END IF
      END DO
      NSOPR=MSOPR

C IPUSED is used later for other purposes, and should be initialized
C to zero.
      DO IPRP=1,NPRPLST
       IPUSED(IPRP)=0
      END DO
C
C PTYPE and SOPRTP is set here if not already set above. Note that this
C is a fallback procedure that should not be used actively. This
C fallback typically assigns the type of WE-reduced density to be
C used for user-specified lists of properties.
C
      DO IPROP=1,NPROP
       IF (PTYPE(IPROP).NE.'UNDEF.  ') CYCLE
       PTYPE(IPROP)='HERMSING'
       IF(PNAME(IPROP).EQ.'VELOCITY') PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'ANGMOM  ') PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'MLTPV  2') PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP)(1:4).EQ.'PSOP') PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'OMQ     ') PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'AMFI    ') PTYPE(IPROP)='ANTITRIP'
       IF(PNAME(IPROP)(1:3).EQ.'ASD') PTYPE(IPROP)='HERMTRIP'
       IF(PNAME(IPROP)(1:6).EQ.'DMP   ') PTYPE(IPROP)='HERMSING'
       IF(PNAME(IPROP).EQ.'EMFR  RA')PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'EMFR  IA')PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'TMOM  RA')PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'TMOM  IA')PTYPE(IPROP)='ANTISING'
       IF(PNAME(IPROP).EQ.'EMFR0  I')PTYPE(IPROP)='ANTITRIP'
       IF(PNAME(IPROP).EQ.'TMOM0  R')PTYPE(IPROP)='HERMTRIP'
       IF(PNAME(IPROP).EQ.'TMOM0  I')PTYPE(IPROP)='ANTITRIP'
       IF(PNAME(IPROP).EQ.'TMOM2  I')PTYPE(IPROP)='ANTISING'
       END DO

      DO ISOPR=1,NSOPR
       IF (SOPRTP(ISOPR).NE.'UNDEF.  ') CYCLE
       SOPRTP(ISOPR)='HERMSING'
       IF(SOPRNM(ISOPR).EQ.'VELOCITY') SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'ANGMOM  ') SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'MLTPV  2') SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR)(1:4).EQ.'PSOP') SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(IPROP).EQ.'OMQ     ') SOPRTP(IPROP)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'AMFI    ') SOPRTP(ISOPR)='ANTITRIP'
       IF(SOPRNM(ISOPR)(1:3).EQ.'ASD') SOPRTP(ISOPR)='HERMTRIP'
       IF(SOPRNM(ISOPR)(1:6).EQ.'DMP   ') SOPRTP(ISOPR)='HERMSING'
       IF(SOPRNM(ISOPR).EQ.'EMFR  RA')SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'EMFR  IA')SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'TMOM  RA')SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'TMOM  IA')SOPRTP(ISOPR)='ANTISING'
       IF(SOPRNM(ISOPR).EQ.'EMFR0  I')SOPRTP(ISOPR)='ANTITRIP'
       IF(SOPRNM(ISOPR).EQ.'TMOM0  R')SOPRTP(ISOPR)='HERMTRIP'
       IF(SOPRNM(ISOPR).EQ.'TMOM0  I')SOPRTP(ISOPR)='ANTITRIP'
       IF(SOPRNM(ISOPR).EQ.'TMOM2  I')SOPRTP(ISOPR)='ANTISING'
      END DO

C Write out various input data:
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
      End Do
*
* determine if there are any matching wavefunctions
      JOBMATCH=.FALSE.
      do job1=1,njob
        do job2=1,job1-1
          if ((MLTPLT(JOB1).EQ.MLTPLT(JOB2)) .AND.
     &        (IRREP(JOB1).EQ.IRREP(JOB2))) then
            JOBMATCH=.TRUE.
          end if
        end do
      end do
* make decision regarding the use of input hamiltonian/diagonal values
      if (ifheff) then
        if (have_heff) then
          DO J=1,NSTATE
            DO I=1,NSTATE
              iadr=(j-1)*nstate+i-1
              iadr2=(i-1)*nstate+j-1
              HAM(i,j)=0.5D0*(Work(L_HEFF+iadr)+
     &                               Work(L_HEFF+iadr2))
            END DO
          END DO
          if (jobmatch) then
            call WarningMessage(1,'HEFF used for a situation where '//
     &        'possible extra interaction between states is ignored!')
          end if
        else
          call WarningMessage(2,'HEFF used but none is available!')
          call Quit_OnUserError
        end if
      else if (ifejob) then
        if (have_heff) then
          call WarningMessage(1,'EJOB used when HEFF is available, '//
     &      'possible extra interaction between states is ignored!')
        end if
        if (have_diag) then
          DO I=1,NSTATE
           HAM(i,i)=Work(LREFENE+i-1)
          END DO
        else if (have_heff) then
          DO I=1,NSTATE
            HAM(i,i)=Work(L_HEFF+(i-1)*nstate+i-1)
          END DO
        else
          call WarningMessage(2,'EJOB used but no energies available!')
          call Quit_OnUserError
        end if
        if (jobmatch) then
          call WarningMessage(1,'EJOB used for a situation where '//
     &      'possible extra interaction between states is ignored!')
        end if
      else if (.not.(ifhext.or.ifhdia.or.ifshft.or.ifhcom)) then
* the user has selected no procedure...
        if (have_heff.and. (.not.jobmatch)) then
          ifheff=.true.
          DO J=1,NSTATE
            DO I=1,NSTATE
              iadr=(j-1)*nstate+i-1
              iadr2=(i-1)*nstate+j-1
              HAM(i,j)=0.5D0*(Work(L_HEFF+iadr)+
     &                               Work(L_HEFF+iadr2))
            END DO
          END DO
        else if (have_diag) then
          ifhdia=.true.
          DO I=0,NSTATE-1
            Work(LHDIAG+I)=Work(LREFENE+i)
          END DO
        end if
      end if
*
      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)'  The following data are common to all the states:'
        WRITE(6,*)'  ------------------------------------------------'
        WRITE(6,*)
     &    '  (note: frozen counts as inactive, deleted as secondary)'
        WRITE(6,*)
        WRITE(6,'(6X,A,I2)')'Nr of irreps:',NSYM
        WRITE(6,*)
        WRITE(6,'(6X,A)')
     &       '           Total     No./Irrep '
        WRITE(6,'(6X,A,8X,8I4)')
     &       'Irrep       ',(I,I=1,NSYM)
        WRITE(6,'(6X,A,8X,8(1X,A))')
     &       '            ',(lIrrep(I),I=1,NSYM)
        WRITE(6,*)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &       'INACTIVE    ',NISHT, (NISH(I),I=1,NSYM)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &       'ACTIVE      ',NASHT, (NASH(I),I=1,NSYM)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &       'SECONDARY   ',NSSHT, (NSSH(I),I=1,NSYM)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &       'BASIS       ',NBST, (NBASF(I),I=1,NSYM)
        WRITE(6,*)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &    'RAS1        ',NRS1T, (NRS1(I),I=1,NSYM)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &    'RAS2        ',NRS2T, (NRS2(I),I=1,NSYM)
        WRITE(6,'(6X,A,I4,4X,8I4)')
     &    'RAS3        ',NRS3T, (NRS3(I),I=1,NSYM)
        WRITE(6,*)
        IF(.NOT.(TRACK.OR.ONLY_OVERLAPS)) THEN
          WRITE(6,*) '      '
     &           //' MATRIX ELEMENTS WILL BE COMPUTED FOR THE FOLLOWING'
     &           //' ONE-ELECTRON OPERATORS, UNLESS ZERO BY SYMMETRY.'
          WRITE(6,*) '  (Herm=Hermitian, Anti=Antihermitian,'//
     &                 ' Sing=Singlet operator, Trip=Triplet operator)'
          Do i1=1,nProp,3
            i2=Min(i1+2,nProp)
            WRITE(6,'(3(5X,A8,1X,I3,1X,A1,A8,A1))')
     &           (PNAME(i),ICOMP(i),'(',PTYPE(i),')',i=i1,i2)
          End Do
        END IF
        IF(IFHAM) THEN
          WRITE(6,*)
          WRITE(6,*)'      EIGENSTATES OF A SPIN-FREE HAMILTONIAN'//
     &            ' WILL BE COMPUTED BASED ON:'
          WRITE(6,*)
* which kind of base hamiltonian is taken?
          IF(IFHEXT) THEN
            WRITE(6,*)' a Hamiltonian matrix that '//
     &                'was supplied in the input.'
          ELSE IF(IFHEFF) THEN
            WRITE(6,*)' a (effective) Hamiltonian matrix that '//
     &                'was read from the wavefunction file(s).'
          ELSE IF(IFEJOB) THEN
            WRITE(6,*)' a Hamiltonian matrix assumed to be diagonal '//
     &          'with energies read from the wavefunction file(s).'
          ELSE
            WRITE(6,*)' A Hamiltonian matrix computed by RASSI.'
          END IF
* which kind of corrections are applied?
          IF(IFHDIA) THEN
            WRITE(6,*)' In addition, the diagonal energies of the '//
     &    'hamiltonian matrix will be replaced by either the user '//
     &    '(HDIAG keyword) or read from the wavefunction file(s).'
          END IF
          IF(IFSHFT) THEN
            WRITE(6,*)' In addition, the diagonal energies of the '//
     &    'hamiltonian matrix will be shifted by the user '//
     &    '(SHIFT keyword).'
          END IF
          IF(NSOPR.GT.0) THEN
            WRITE(6,*)' SO coupling elements will be added.'
            WRITE(6,*)'      EIGENSTATES OF SPIN-ORBIT HAMILTONIAN'//
     &            ' WILL BE COMPUTED'
          END IF
        END IF
        IF(.NOT.(TRACK.OR.ONLY_OVERLAPS)) THEN
          IF(NSOPR.GT.0) THEN
          WRITE(6,*) '       MATRIX ELEMENTS OVER SPIN EIGENSTATES FOR:'
          Do i1=1,NSOPR,3
            i2=Min(i1+2,NSOPR)
            WRITE(6,'(3(5X,A8,1X,I3,1X,A1,A8,A1))')
     &           (SOPRNM(i),ISOCMP(i),'(',SOPRTP(i),')',i=i1,i2)
          End Do
          END IF
        END IF
      END IF
      IF(IPGLOB.GE.DEBUG) THEN
        WRITE(6,*)'Initial default flags are:'
        WRITE(6,*)'     PRSXY :',PRSXY
        WRITE(6,*)'     PRORB :',PRORB
        WRITE(6,*)'     PRTRA :',PRTRA
        WRITE(6,*)'     PRCI  :',PRCI
        WRITE(6,*)'     IFHAM :',IFHAM
        WRITE(6,*)'     IFHEXT:',IFHEXT
        WRITE(6,*)'     IFHEFF:',IFHEFF
        WRITE(6,*)'     IFEJOB:',IFEJOB
        WRITE(6,*)'     IFHCOM:',IFHCOM
        WRITE(6,*)'     IFSHFT:',IFSHFT
        WRITE(6,*)'     IFHDIA:',IFHDIA
        WRITE(6,*)'     IFSO  :',IFSO
        WRITE(6,*)'     NATO  :',NATO
        WRITE(6,*)'     RFPERT:',RFPERT
        WRITE(6,*)'     TOFILE:',ToFile
        WRITE(6,*)'     PRXVR :',PRXVR
        WRITE(6,*)'     PRXVE :',PRXVE
        WRITE(6,*)'     PRXVS :',PRXVS
        WRITE(6,*)'     PRMER :',PRMER
        WRITE(6,*)'     PRMEE :',PRMEE
        WRITE(6,*)'     PRMES :',PRMES
      END IF
      IF(IPGLOB.GE.USUAL) THEN
       IF(NATO.AND.(NRNATO.GT.0)) THEN
        WRITE(6,*)' Natural orbitals will be computed for the'
        WRITE(6,*)' lowest eigenstates. NRNATO=',NRNATO
       END IF
       IF(BINA) THEN
        WRITE(6,*)' Bi-natural orbitals will be computed for the'
        WRITE(6,*)' following pairs of states:'
        WRITE(6,'(5X,8(2X,A1,I2,A1,I2,A1))')
     &             ('(',IBINA(1,I),',',IBINA(2,I),')',I=1,NBINA)
       END IF
       WRITE(6,*)
       WRITE(6,*)' Nr of states:',NSTATE
       DO II=1,NSTATE,20
        III=MIN(II+19,NSTATE)
        WRITE(6,*)
        WRITE(6,'(1X,A8,5x,20I4)')'  State:',(I,I=II,III)
        WRITE(6,'(1X,A8,5x,20I4)')' JobIph:',
     &                             (iWork(lJBNUM+I-1),I=II,III)
        WRITE(6,'(1X,A8,5x,20I4)')'Root nr:',
     &                             (iWork(lLROOT+I-1),I=II,III)
       END DO
       IF(IFSHFT) THEN
         WRITE(6,*)
         WRITE(6,*)'Each input state will be shifted with an individual'
         WRITE(6,*)'amount of energy. These energy shifts are (a.u.):'
         WRITE(6,'(1X,5F16.8)')(ESHFT(I),I=1,NSTATE)
       END IF
      END IF

C Added by Ungur Liviu on 04.11.2009
C Addition of NSTATE, JBNUM, and LROOT to RunFile.

       CALL Put_iscalar('NSTATE_SINGLE',NSTATE)
       CALL Put_iArray('JBNUM_SINGLE',iWork(lJBNUM),NSTATE)
       CALL Put_iArray('LROOT_SINGLE',iWork(lLROOT),NSTATE)
*
* Generate the quadrature points for isotropic integration of the exponential operator
*
      If (Do_TMOM) Then
        If (Do_SK) Then
          nQuad=1
        Else
          nk_Vector = 1
          Call Setup_O()
          Call Do_Lebedev_Sym(L_Eff,nQuad,ipR)
          Call Free_O()
          Call Free_Work(ipR)
        End If
      Else
        nk_Vector = 0
        nQuad = 0
      End If
*
      CALL XFLUSH(6)
      CALL QEXIT(ROUTINE)
      RETURN
      END
