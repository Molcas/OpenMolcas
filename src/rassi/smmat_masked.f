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
      SUBROUTINE SMMAT_MASKED(PROP,PRMAT,NSS,PRLBL,IPRCMP,ISS_INDEX,
     &                        IST,INUM,JST,JNUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) PRLBL
      DIMENSION PRMAT(NSS,NSS)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SMMAT_MASKED')
      Integer INUM, JNUM
      Integer ISS_INDEX(NSTATE+1), IST(INUM), JST(JNUM)
#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP)
      REAL*8, EXTERNAL :: DCLEBS
*
      IPRNUM=-1
C IFSPIN takes values the values 0,1,2
C 0 = spin free property
C 1 = spin operator (S)
C 2 = spin dependent property, triplet operator
      IFSPIN=0

      DO IPROP=1,NPROP
         IF (PRLBL.EQ.PNAME(IPROP)) THEN
            IF (PRLBL(1:5).eq.'TMOM0') THEN
               IFSPIN=2
               IPRNUM=IPROP
               EXIT
            ELSE
               IFSPIN=0
               IF (IPRCMP.EQ.ICOMP(IPROP)) THEN
                  IPRNUM=IPROP
                  EXIT
               END IF
            END IF
         ELSE IF (PRLBL(1:4).EQ.'SPIN') THEN
            IFSPIN=1
            IPRNUM=0
            EXIT
         END IF
      END DO
      IF (IPRNUM.EQ.-1) THEN
         Write (6,*) 'SMMAT_MASKED, Abend IPRNUM.EQ.-1'
         Write (6,*) 'SMMAT_MASKED, PRLBL=','>',PRLBL,'<'
         Call Abend()
      ENDIF

C Mapping from spin states to spin-free state and to spin:
      DO I=1,INUM
         ISTATE=IST(I)
         ISS=ISS_INDEX(ISTATE)
         MPLET1=ISS_INDEX(ISTATE+1)-ISS_INDEX(ISTATE)
         S1=0.5D0*(MPLET1-1)
         DO MSPROJ1=-MPLET1+1,MPLET1-1,2
            SM1=0.5D0*MSPROJ1
            ISS=ISS+1

            DO J=1,JNUM
               JSTATE=JST(J)
               JSS=ISS_INDEX(JSTATE)
               MPLET2=ISS_INDEX(JSTATE+1)-ISS_INDEX(JSTATE)
               S2=0.5D0*(MPLET2-1)
               DO MSPROJ2=-MPLET2+1,MPLET2-1,2
                  SM2=0.5D0*MSPROJ2
                  JSS=JSS+1

                  IF (IFSPIN.EQ.0 .AND. IPRNUM.NE.0) THEN
                     IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
                        PRMAT(ISS,JSS)=PROP(ISTATE,JSTATE,IPRNUM)
                     ELSE
                        PRMAT(ISS,JSS)=0.0D0
                     END IF
                  ELSE IF (IFSPIN.EQ.1 .AND. IPRNUM.EQ.0) THEN
                     SXMER=0.0D0
                     SYMEI=0.0D0
                     SZMER=0.0D0
                     SMINUS=0.0D0
                     SPLUS=0.0D0
                     IF((ISTATE.EQ.JSTATE).AND.(MPLET1.eq.MPLET2)) THEN
                        IF (MSPROJ1.eq.MSPROJ2-2) THEN
                           SMINUS=SQRT((S1+SM2)*(S1-SM1))
                           SXMER= 0.5D0*SMINUS
                           SYMEI= 0.5D0*SMINUS
                        ELSE IF (MSPROJ1.eq.MSPROJ2) THEN
                           SZMER=SM1
                        ELSE IF (MSPROJ1.eq.MSPROJ2+2) THEN
                           SPLUS=SQRT((S1+SM1)*(S1-SM2))
                           SXMER= 0.5D0*SPLUS
                           SYMEI=-0.5D0*SPLUS
                        END IF
                        IF (IPRCMP.EQ.1) THEN
                           PRMAT(ISS,JSS)=SXMER
                        ELSE IF (IPRCMP.EQ.2) THEN
                           PRMAT(ISS,JSS)=SYMEI
                        ELSE IF (IPRCMP.EQ.3) THEN
                           PRMAT(ISS,JSS)=SZMER
                        END IF
                     ELSE
                        PRMAT(ISS,JSS)=0.0D0
                     END IF
                  ELSE IF (IFSPIN.EQ.2) THEN
C
C                 The code here is a replica from smmat.f. Look in
C                 that source for comments.
C
                     FACT=1.0D0/SQRT(DBLE(MPLET1))
                     IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
                     CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
                     CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
                     CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
                     CGX= SQRT(0.5D0)*(CGM-CGP)
                     CGY=-SQRT(0.5D0)*(CGM+CGP)
*
                     EXPKR=PROP(ISTATE,JSTATE,IPRNUM)
*
                     IF (IPRCMP.EQ.1) THEN
                        EXPKR=EXPKR*CGX
                     ELSE IF (IPRCMP.EQ.2) THEN
                        EXPKR=EXPKR*CGY
                     ELSE IF (IPRCMP.EQ.3) THEN
                        EXPKR=EXPKR*CG0
                     END IF
                     PRMAT(ISS,JSS)= EXPKR
                  END IF

               END DO !MSPROJ2
            END DO !JSTATE
         END DO !MSPROJ1
      END DO !ISTATE

      RETURN
      END
