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
      SUBROUTINE SMMAT(PROP,PRMAT,NSS,PRLBL,IPRCMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) PRLBL
      DIMENSION PRMAT(NSS,NSS)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SMMAT')
#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP)
*
      IPRNUM=0
C IFSPIN takes values the values 0,1,2
C 0 = spin free property
C 1 = spin operator (S)
C 2 = spin dependent property, triplet operator
      IFSPIN=0

      DO IPROP=1,NPROP
        IF (PRLBL.EQ.PNAME(IPROP)) THEN
           IFSPIN=0
           IF(IPRCMP.EQ.ICOMP(IPROP)) IPRNUM=IPROP
        ELSE IF (PRLBL(1:4).EQ.'SPIN') THEN
           IFSPIN=1
        ELSE IF (PRLBL(1:5).EQ.'TMOM0') THEN
           IFSPIN=2
*
*          Note that the integral is complex. Select the real or the
*          imaginary component here.
*
           IF (PRLBL(1:8).EQ.PNAME(IPROP)) IPRNUM=IPROP
        END IF
      END DO
C Mapping from spin states to spin-free state and to spin:
      ISS=0
      DO ISTATE=1,NSTATE
       JOB1=iWork(lJBNUM+ISTATE-1)
       MPLET1=MLTPLT(JOB1)
       S1=0.5D0*DBLE(MPLET1-1)

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
        SM1=0.5D0*DBLE(MSPROJ1)
        ISS=ISS+1

        JSS=0

        DO JSTATE=1,NSTATE
         JOB2=iWork(lJBNUM+JSTATE-1)
         MPLET2=MLTPLT(JOB2)
         S2=0.5D0*DBLE(MPLET2-1)

         DO MSPROJ2=-MPLET2+1,MPLET2-1,2
          SM2=0.5D0*DBLE(MSPROJ2)
          JSS=JSS+1

          IF (IFSPIN.EQ.0 .AND. IPRNUM.NE.0) THEN
                  IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
                          PRMAT(ISS,JSS)=PROP(ISTATE,JSTATE,IPRNUM)
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
                  END IF
          ELSE IF (IFSPIN.EQ.2) THEN
C 1-electron triplet operator so only Delta S =0,+-1 and Delta MS =0,+-1
C Notice S1,S2 and SM1,SM2 need not to be integers
C Hence MPLET1,MPLET2 and MSPROJ1,MSPROJ2 are used
C Notice SMINUS and SPLUS is interchanged compared to above
C Notice that the Y part is imaginary
C
C see section 3 (Spin_orbit coupling in RASSI) in
C P A Malmqvist, et. al CPL, 357 (2002) 230-240
C for details
C
C Note that we work on the x, y, and z components at this time.
C
C On page 234 we have the notation V^{AB}(x), that is the
C potential has Cartesian components. Here, however, this is
C partitioned in a slightly different way since we have that
C V^{AB}(x)=(k x e_l)_x V^{AB}. We will only handle the
C V^{AB} part.
C
C Hence, we will compute the contributions to T(i), i=x,y,z
C here and form the inner product
C (k x e_l)_i V^{AB} . T(i)
C outside the code.
C
                  IF(ABS(MPLET1-MPLET2).GT.2) CYCLE
                  IF(ABS(MSPROJ1-MSPROJ2).GT.2) CYCLE
C
                  SXMER =0.0D0
                  SYMEI =0.0D0
                  SZMER =0.0D0
                  SMINUS=0.0D0
                  SPLUS =0.0D0
                  ONE   =1.0D0
                  TWO   =2.0D0
C
                  IF(MPLET1+2.EQ.MPLET2) THEN ! <SM|O|S+1M+?>
C
C                   MSPROJ1-MSPROJ2=-2
                    IF (MSPROJ1+2.eq.MSPROJ2) THEN ! <SM|O|S+1M+1>
                      SPLUS =-0.5D0*SQRT((S1+SM1+ONE)*(S1+SM1+TWO))
                      SXMER =+SPLUS
                      SYMEI =+SPLUS
C
C                   MSPROJ1-MSPROJ2= 0
                    ELSE IF (MSPROJ1.eq.MSPROJ2) THEN ! <SM|O|S+1M>
                      SZMER =SQRT((S1+ONE)**2-SM1**2)
C
C                   MSPROJ1-MSPROJ2=+2
                    ELSE IF (MSPROJ1-2.eq.MSPROJ2) THEN ! <SM|O|S+1M-1>
                      SMINUS=-0.5D0*SQRT((S1-SM1+ONE)*(S1-SM1+TWO))
                      SXMER =-SMINUS
                      SYMEI =+SMINUS
                    END IF
C
                  ELSE IF(MPLET1.EQ.MPLET2) THEN ! <SM|O|SM+?>
C
C                   MSPROJ1-MSPROJ2=-2
                    IF (MSPROJ1+2.eq.MSPROJ2) THEN ! <SM|O|SM+1>
                      SPLUS = 0.5D0*SQRT((S1-SM1)*(S1+SM1+ONE))
                      SXMER =+SPLUS
                      SYMEI =+SPLUS
C
C                   MSPROJ1-MSPROJ2= 0
                    ELSE IF (MSPROJ1.eq.MSPROJ2) THEN ! <SM|O|SM>
                      SZMER =SM1
C
C                   MSPROJ1-MSPROJ2=+2
                    ELSE IF (MSPROJ1-2.eq.MSPROJ2) THEN ! <SM|O|SM-1>
                      SMINUS=-0.5D0*SQRT((S1+SM1)*(S1-SM1+ONE))
                      SXMER =-SMINUS
                      SYMEI =+SMINUS
                    END IF
C
                  ELSE IF(MPLET1-2.EQ.MPLET2) THEN ! <SM|O|S-1M+?>
C
C                   MSPROJ1-MSPROJ2=-2
                    IF (MSPROJ1+2.eq.MSPROJ2) THEN ! <SM|O|S-1M+1>
                      SPLUS =0.5D0*SQRT((S1-SM1)*(S1-SM1-ONE))
                      SXMER = SPLUS
                      SYMEI = SPLUS
C
C                   MSPROJ1-MSPROJ2= 0
                    ELSE IF (MSPROJ1.eq.MSPROJ2) THEN ! <SM|O|S-1M>
                      SZMER =SQRT(S1**2-SM1**2)
C
C                   MSPROJ1-MSPROJ2=+2
                    ELSE IF (MSPROJ1-2.eq.MSPROJ2) THEN ! <SM|O|S-1M-1>
                      SMINUS=0.5D0*SQRT((S1+SM1)*(S1+SM1-ONE))
                      SXMER =-SMINUS
                      SYMEI = SMINUS
                    END IF
C
                  END IF
C
                  IF (IPRCMP.EQ.1) THEN
                    PRMAT(ISS,JSS)=SXMER * PROP(ISTATE,JSTATE,IPRNUM)
                  ELSE IF (IPRCMP.EQ.2) THEN
                    PRMAT(ISS,JSS)=SYMEI * PROP(ISTATE,JSTATE,IPRNUM)
                  ELSE IF (IPRCMP.EQ.3) THEN
                    PRMAT(ISS,JSS)=SZMER * PROP(ISTATE,JSTATE,IPRNUM)
                  END IF
          END IF

         END DO
        END DO
       END DO
      END DO

      RETURN
      END
