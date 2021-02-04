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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE SIZES()
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"


C Available workspace right now:
      CALL GETMEM('LDUM','MAX','REAL',LDUM,MXLEFT)
C 250000 words margin, for various purposes.
      NBOTTOM=250000
C POLY package:
C MINBUF=Min acceptable nr of buffers in POLY3 for efficiency.
      MINBUF=1
      NPOLY=NBOTTOM
C Sizes of G10,G20,G02,G12,G30: Faked, since no more used.
      NG10=1
      NG20=1
      NG02=1
      NG30=1
      NG12=1

      IF(NACTEL.GT.0) THEN
        NGARR=2*(NG10+NG20+NG02+NG30+NG12)
        NPOLY=NPOLY+NGARR+5*MXCI
        MXLFT=MXLEFT-NPOLY
        XX=0.5D0*DBLE(3*MXCI+3)
        YY=DBLE(MXLFT)
C Allowed by memory:
        NPLBUF=INT(SQRT(YY+XX**2)-XX)
C Preferred size for efficiency:
        NPLBUF=MAX(MINBUF,NPLBUF)
C Actual max needed in-core:
        NPLBUF=MIN(NPLBUF,(NASHT*(NASHT+1))/2)
        NPOLY=NPOLY+NPLBUF*(NPLBUF+3+3*MXCI)
*        WRITE(*,*)' Memory requirements for POLY3 (Above SGUGA):'
*        WRITE(*,*)'   CI vector                      :',MXCI
*        WRITE(*,*)'   Extra margin for small scratch :',NBOTTOM
*        WRITE(*,*)'   Ten arrays G10,G20..F12        :',NGARR
*        WRITE(*,*)'   H0,SGM0,C1,C2                  :',4*MXCI
*        WRITE(*,*)'   VXYZB,TUBUF,DTUB               :',3*NPLBUF
*        WRITE(*,*)'   A,B1,B2                        :',3*NPLBUF*MXCI
*        WRITE(*,*)'   SCR                            :',NPLBUF**2
*        WRITE(*,*)
*        WRITE(*,*)' Total, for POLY3                 :',NPOLY
*        WRITE(*,*)
      END IF

C Precompute sizes, offsets etc.
      CALL SUPINI

C SBMAT need:
*     NG3C=0
*     DO ISYM=1,NSYM
*       N=NTUV(ISYM)
*       NG3C=NG3C+(N*(N+1))/2
*     END DO
*     NG3C=iPARDIV(NG3TOT,NG2)

C Sizes and addresses to lists:
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        ISL2=MUL(ISL1,ISL3)
        NLIST(ISL1,ISL3,1)= NASH(ISL2)*NTU(ISL3)*2
        NLIST(ISL1,ISL3,2)= NLIST(ISL1,ISL3,1)
        NLIST(ISL1,ISL3,3)= NASH(ISL2)*NTU(ISL3)
        NLIST(ISL1,ISL3,4)= NASH(ISL2)*NTGTU(ISL3)*2
        NLIST(ISL1,ISL3,5)= NLIST(ISL1,ISL3,3)
        NLIST(ISL1,ISL3,6)= NLIST(ISL1,ISL3,4)
        NLIST(ISL1,ISL3,7)= NASH(ISL2)*NASH(ISL3)*2
        NLIST(ISL1,ISL3,8)= NLIST(ISL1,ISL3,7)
        NLIST(ISL1,ISL3,9)= NASH(ISL2)*NASH(ISL3)
        NLIST(ISL1,ISL3,10)= NLIST(ISL1,ISL3,9)
        IF(ISL1.EQ.1) NLIST(ISL1,ISL3,10)
     &   =NLIST(ISL1,ISL3,9)-NASH(ISL2)
        NLIST(ISL1,ISL3,12)= NASH(ISL1)*NASH(ISL2)
        NLIST(ISL1,ISL3,13)= NLIST(ISL1,ISL3,12)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,13)
     &   =NLIST(ISL1,ISL3,12)-NASH(ISL1)
        NLIST(ISL1,ISL3,11)=NLIST(ISL1,ISL3,12)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,11)
     &   =NLIST(ISL1,ISL3,11)+NASH(ISL1)*NASHT
        NLIST(ISL1,ISL3,14)= NISH(ISL1)*NISH(ISL2)
        NLIST(ISL1,ISL3,15)= NLIST(ISL1,ISL3,14)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,15)
     &   =NLIST(ISL1,ISL3,14)-NISH(ISL1)
        NLIST(ISL1,ISL3,16)= NSSH(ISL1)*NSSH(ISL2)
        NLIST(ISL1,ISL3,17)= NLIST(ISL1,ISL3,16)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,17)
     &   =NLIST(ISL1,ISL3,16)-NSSH(ISL1)
       END DO
      END  DO
      NLISTS=0
      DO ILIST=1,17
       DO ISL1=1,NSYM
        DO ISL3=1,NSYM
         NLISTS=NLISTS+NLIST(ISL1,ISL3,ILIST)
        END DO
       END DO
      END DO
      NLSTOT=4*NLISTS

C maximum orbitals in an irrep
      NOMAX=0
      DO ISYM=1,NSYM
        NOMAX=MAX(NOMAX,NORB(ISYM))
      END DO
      IF (.NOT.IFCHOL) THEN
C MKRHS needs:
        NMKRHS=2*NOMAX**2
        NMX=0
        DO ICASE=1,13
          DO ISYM=1,NSYM
            NIN=NINDEP(ISYM,ICASE)
            IF(NIN.EQ.0) GOTO 10
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NV=NAS*NIS
            IF(NV.EQ.0) GOTO 10
            IF(ICASE.GT.11) THEN
              N=NV
            ELSE
              N=2*NV+(NAS*(NAS+1))/2
            END IF
            NMX=MAX(N,NMX)
  10        CONTINUE
          END DO
        END DO
        NMKRHS=NMKRHS+NMX
*       IF(MAXIT.EQ.0) THEN
*         NSIGMA=0
*       ELSE
*         NSIGMA=NBOTTOM+NLSTOT+MMX
*       END IF
      ELSE
C RHSALL2 and ADDRHS needs:
        NMKRHS=0
        DO ICASE=1,13
          DO ISYM=1,NSYM
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NMKRHS=MAX(NMKRHS,NAS*NIS)
          END DO
        END DO
        NMKRHS = iPARDIV(NMKRHS,2*NOMAX**2)
        NMKRHS = 2*NMKRHS
      END IF
      NMKRHS=NMKRHS+NBOTTOM

C PCG/(new)SIGMA routine needs: twice a global array RHS size for
C transformations (overrides NSIGMA computed above)
      NSIGMA_OUTER=0
      NSIGMA_INNER=0
      NVCUTIL=0
      DO ICASE=1,13
      DO ISYM=1,NSYM
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NRHSP = iPARDIV(NAS*NIS,0)
      IF (ICASE.GT.11) THEN
        NSIGMA_INNER=MAX(NSIGMA_INNER,NRHSP)
      ELSE
        NSIGMA_INNER=MAX(NSIGMA_INNER,NAS*NIS+NRHSP)
        NSIGMA_OUTER=MAX(NSIGMA_OUTER,NAS*NIS+NRHSP)
      END IF
      NVCUTIL=MAX(NVCUTIL,2*NRHSP)
      END DO
      END DO
      NSIGMA=MAX(NSIGMA_INNER+NSIGMA_OUTER,NVCUTIL)
      NSIGMA=NSIGMA+NBOTTOM

C PRPCTL needs:
C In DIADNS alone, NDD words are needed:
#ifdef _DEBUGPRINT_
      WRITE(6,*)' Memory requirements for PRPCTL (Above SGUGA).'
      WRITE(6,*)
      WRITE(6,*)' PRP1) First phase of PRPCTL.'
      WRITE(6,*)
      WRITE(6,*)'    A) DIADNS alone, broken up in case/symm:'
#endif
      MMX=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 11
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NV=NAS*NIS
          IF(NV.EQ.0) GOTO 11
          IF(ICASE.GT.11) THEN
            M=2*NV
          ELSE
            M=3*NV+(NAS*(NAS+1))/2
          END IF
          MMX=MAX(M,MMX)
  11      CONTINUE
        END DO
      END DO

      NDD=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NX=0
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            NIS=NISUP(ISYM,ICASE)
            IF(NIS.GT.0) THEN
              NAS=NASUP(ISYM,ICASE)
#ifdef _DEBUGPRINT_
              WRITE(6,*)' Case, Symm:',ICASE,ISYM
              WRITE(6,*)' NIN,NAS,NIS:',NIN,NAS,NIS
              WRITE(6,*)' NIMX,NSMX:',NIMX,NAMX
#endif
              IF(ICASE.EQ.2 .OR. ICASE.EQ.3) NX=NIN*NIMX**2
              IF(ICASE.EQ.6 .OR. ICASE.EQ.7) NX=NIN*NSMX*NIMX**2
              IF(ICASE.EQ.8 .OR. ICASE.EQ.9) NX=NIN*NSMX**2
              IF(ICASE.EQ.10 .OR. ICASE.EQ.11) NX=NIN*NIMX*NSMX**2
              IF(ICASE.EQ.12 .OR. ICASE.EQ.13) THEN
                NX=MAX(NAS*NIMX**2,NIS*NSMX**2)
              END IF
            END IF
          END IF
#ifdef _DEBUGPRINT_
      WRITE(6,'(1x,a,2i4,5x,i8)')'      Case, symm:',ICASE,ISYM,2*NX
#endif
          NDD=MAX(2*NX,NDD)
        END DO
      END DO

#ifdef _DEBUGPRINT_
      WRITE(6,*)'       DIADNS needs the maximum, or NDD=',NDD
#endif
      NCMO=NBSQT
      NPRP1=NBOTTOM+NCMO+notri+NLSTOT+2*NOSQT+MMX+NDD

#ifdef _DEBUGPRINT_
      WRITE(6,*)
      WRITE(6,*)'    B) Also needed, for 1st phase of PRPCTL:'
      WRITE(6,'(1x,a,i8)')'       NBOTTOM:',NBOTTOM
      WRITE(6,'(1x,a,i8)')'       NCMO   :',NCMO
      WRITE(6,'(1x,a,i8)')'       notri :',notri
      WRITE(6,'(1x,a,i8)')'       NLSTOT :',NLSTOT
      WRITE(6,'(1x,a,i8)')'     2*NOSQT  :',2*NOSQT
      WRITE(6,'(1x,a,i8)')'       MMX    :',MMX
#endif
      NPRP2=NBOTTOM+2*NCMO+notri+NBAST*(NBAST+1)

#ifdef _DEBUGPRINT_
      WRITE(6,*)' PRP2) Second phase of PRPCTL.'
      WRITE(6,*)
      WRITE(6,'(1x,a,i8)')'       NBOTTOM:',NBOTTOM
      WRITE(6,'(1x,a,i8)')'     2*NCMO   :',2*NCMO
      WRITE(6,'(1x,a,i8)')'       notri :',notri
      WRITE(6,'(1x,a,i8)')'NBAST*(NBAST+1)',NBAST*(NBAST+1)
#endif

      NPRP=0
      IF (IFPROP) NPRP=MAX(NPRP1,NPRP2)
      IF ( IPRGLB.GE.USUAL) THEN
        WRITE(6,'(20A4)')('----',I=1,20)
        WRITE(6,*)'Estimated memory requirements:'
        WRITE(6,'(a,i12)')'  POLY3 :             ',NPOLY
        WRITE(6,'(a,i12)')'  RHS:                ',NMKRHS
        WRITE(6,'(a,i12)')'  SIGMA :             ',NSIGMA
        WRITE(6,'(a,i12)')'  PRPCTL:             ',NPRP
*SVC: NPRP includes NDD if it is needed, so this is confusing
*       WRITE(6,'(a,i12)')'  DIADNS:             ',NDD
        WRITE(6,'(a,i12)')' Available workspace: ',MXLEFT
        WRITE(6,*)
      ENDIF

      NEED0=MAX(NPOLY,NMKRHS,NSIGMA)
      NEED=MAX(NEED0,NPRP)
      IF(NEED.GT.MXLEFT) THEN
       IF(NEED0.LE.MXLEFT) THEN
        WRITE(6,'(5X,26A)') ('*',i=1,26)
        WRITE(6,'(5X,A)')' Memory problem!! The memory is insufficient '
        WRITE(6,'(5X,A)')' for the property section.'
        WRITE(6,'(5X,A,I5,A)')'* Need at least ',2+NEED/119000,' MB *'
        WRITE(6,'(5X,A)')' The property section will be skipped.'
        WRITE(6,'(5X,26A)') ('*',i=1,26)
        NEED=NEED0
        IFPROP=.False.
       ELSE
        IF (.NOT.IFCHOL) THEN
C not a Cholesky calculation, keep old memory requirements
         WRITE(6,'(5X,26A)') ('*',i=1,26)
         WRITE(6,'(5X,A)')'* Insufficient memory !! *'
         WRITE(6,'(5X,A,I5,A)')'* Need at least ',2+NEED0/119000,' MB *'
         WRITE(6,'(5X,26A)') ('*',i=1,26)
         WRITE(6,*)
         WRITE(6,*)' Program execution stops -- sorry!'
         WRITE(6,*)
         CALL Quit(_RC_MEMORY_ERROR_)
        ELSE
C This is a Cholesky calculation, only give recommended amount
         IF (MAX(NPOLY,NSIGMA).LE.MXLEFT) THEN
          WRITE(6,'(5X,40A)') ('*',i=1,40)
          WRITE(6,'(5X,A,I5,A)')'* Below comfortable memory of ',
     &                           2+NEED0/119000,' MB *'
          WRITE(6,'(5X,40A)') ('*',i=1,40)
          WRITE(6,*)
          WRITE(6,*)' Program will try to adapt'
          WRITE(6,*)' If it fails, please raise the memory to at least',
     &            2+MAX(NPOLY,INT(0.55D0*NMKRHS),NSIGMA)/119000,' MB'
          WRITE(6,*)' (Maybe more, this aint rocket science)'
          WRITE(6,*)
          WRITE(6,*)' With print level DEBUG you will get some more'
          WRITE(6,*)' informative memory estimates during computation'
          WRITE(6,*)
         ELSE
          WRITE(6,'(5X,26A)') ('*',i=1,26)
          WRITE(6,'(5X,A)')'* Insufficient memory... *'
          WRITE(6,'(5X,26A)') ('*',i=1,26)
          WRITE(6,*)
          WRITE(6,'(A,I6,A)')' If possible, I would like to have   ',
     &                NEED0/119000, ' MB'
          WRITE(6,'(A,I6,A)')' Please raise the memory to at least ',
     &            MAX(NPOLY,INT(0.55D0*NMKRHS),NSIGMA)/119000,' MB'
          WRITE(6,*)' (Maybe more, this aint rocket science)'
          WRITE(6,*)
          CALL Quit(_RC_MEMORY_ERROR_)
         END IF
        END IF
       END IF
      END IF


      RETURN
      END
