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
* Copyright (C) 2015, Roland Lindh                                     *
************************************************************************
      SUBROUTINE MK_PROP(PROP,IPROP,ISTATE_,JSTATE_,LABEL,ITYPE,
     &                   BUFF,NBUFF,DENS,NDENS,MASK,ISY12,IOFF)
      use OneDat, only: sOpSiz
      use Constants, only: Zero, Two
      use Cntrl, only: NSTATE, NPROP, PNUC, PORIG, ICOMP, IPUSED, PNAME
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      IMPLICIT None
************************************************************************
*     Objective: to compute the transition property between state      *
*                ISTATE and JSTATE of property IPROP.                  *
*                                                                      *
*     This routine will be generalized to a direct routine later.      *
*                                                                      *
*     Author: Roland Lindh, Uppsala University, 23 Dec. 2015           *
************************************************************************
#include "rassi.fh"
      INTEGER IPROP, ISTATE_, JSTATE_, ITYPE, NBUFF, NDENS, MASK, ISY12
      CHARACTER(LEN=8) LABEL
      REAL*8 PROP(NSTATE,NSTATE,NPROP), BUFF(NBUFF), DENS(NDENS,4)
      INTEGER IOFF(8)

      Integer IDUM(1), ISTATE, JSTATE, IC, IOPT, NSIZ, IRC, ISCHK, IINT,
     &        ISY1, NB1, ISY2, NB2, I12, NB12, IPOS
      REAL*8 PSUM
      REAL*8, External:: DDot_

      ISTATE=Max(ISTATE_,JSTATE_)
      JSTATE=Min(ISTATE_,JSTATE_)
*
      IF (LABEL(1:4).eq.'ASD ') LABEL(1:5)='MAGXP'
      IC=ICOMP(IPROP)
!     Write (*,*) 'Mk_Prop: Label=',Label
!     Write (*,*) 'Mk_Prop:    IC=',IC
      IOPT=ibset(0,sOpSiz)
      NSIZ=0
      CALL iRDONE(IRC,IOPT,LABEL,IC,IDUM,ISCHK)
      IF(IRC.eq.0) NSIZ=IDUM(1)
      IF(MOD(ISCHK/MASK,2).EQ.0) RETURN
      IOPT=0
C Rulin: The 'spin-dependent' part of hyperfine contribution
      IF (LABEL(1:5).eq.'MAGXP') THEN
        CALL HFCSD(LABEL,IC,BUFF,NBUFF,NSIZ,ISCHK)
        LABEL(1:5) = 'ASD  '
      ELSE
        CALL RDONE(IRC,IOPT,LABEL,IC,BUFF,ISCHK)
      END IF
!     Write (*,*) 'NBUFF,NSIZ=',NBUFF,NSIZ
      IF ( IRC.NE.0.AND.LABEL(1:4).NE.'TMOM' ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)')'*** ERROR IN SUBROUTINE MK_PROP ***'
         WRITE(6,'(6X,A)')'  FAILED IN READING FROM  ONEINT'
         WRITE(6,'(6X,A,A)')'  LABEL     = ',LABEL
         WRITE(6,'(6X,A,I2)')'  COMPONENT = ',IC
         WRITE(6,*)
         RETURN
      END IF
      IPUSED(IPROP)=1
C IF THIS IS THE FIRST CALL TO THE SUBROUTINE, PICK UP SOME DATA:
*       IF(ICALL.EQ.0) THEN
C-SVC: for safety reasons, always pick up the data, since it is not
C      necessarily done on the first call!
C PICK UP THE ORIGIN COORDINATES:
      PORIG(1,IPROP)=BUFF(NSIZ+1)
      PORIG(2,IPROP)=BUFF(NSIZ+2)
      PORIG(3,IPROP)=BUFF(NSIZ+3)
C PICK UP THE NUCLEAR CONTRIBUTION FROM INTEGRAL BUFFER
      IF (PNAME(IPROP)(1:3).EQ.'ASD'.OR.
     &    PNAME(IPROP)(1:3).EQ.'PSO') THEN
         Write(6,*) "Removing nuclear contrib from ASD and PSO: "
         PNUC(IPROP)=0.0d0
      ELSE IF (ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
         PNUC(IPROP)=0.0d0
      ELSE
         PNUC(IPROP)=BUFF(NSIZ+4)
      END IF
      IINT=1
      PSUM=Zero
      DO ISY1=1,NSYM
         NB1=NBASF(ISY1)
         IF(NB1.EQ.0) cycle
         DO ISY2=1,ISY1
            I12=MUL(ISY1,ISY2)
            IF(IAND(2**(I12-1),ISCHK).EQ.0) cycle
            NB2=NBASF(ISY2)
            IF(NB2.EQ.0) cycle
            NB12=NB1*NB2
            IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
            IF (I12.EQ.ISY12) THEN
               IPOS=IOFF(ISY1)+1
               PSUM=PSUM+DDOT_(NB12,BUFF(IINT),1,DENS(IPOS,ITYPE),1)
            END IF
            IINT=IINT+NB12
        END DO
      END DO
C     Write (*,*) 'PSUM=',PSUM,LABEL,IC
C IN THE CASE OF MULTIPOLES, CHANGE SIGN TO ACCOUNT FOR THE NEGATIVE
C ELECTRONIC CHARGE AS COMPARED TO THE NUCLEAR CONTRIBUTION.
      IF(LABEL(1:5).EQ.'MLTPL') PSUM=-PSUM
C In the case of AMFI integrals, they should be multiplied by 2.
C The reason for the factor of two is that this program uses spin
C (uses Clebsch-Gordan coefficients to define Wigner-Eckart
C  reduced matrix elements of spin-tensor properties)
C while the AMFI authors used Pauli matrices.
      IF(LABEL(1:4).EQ.'AMFI') PSUM=Two*PSUM
      PROP(ISTATE,JSTATE,IPROP)=PSUM
      IF (ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
         PROP(JSTATE,ISTATE,IPROP)=PSUM
      ELSE
         PROP(JSTATE,ISTATE,IPROP)=-PSUM
      END IF
*
      END SUBROUTINE MK_PROP
