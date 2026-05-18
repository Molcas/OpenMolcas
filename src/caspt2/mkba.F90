!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
!*******************************************************************************
! Case A (ICASE=1)
!*******************************************************************************
      SUBROUTINE MKBA(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)
      use definitions, only: iwp, wp, Byte
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM, NINDEP, NTUV
      IMPLICIT NONE
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER(KIND=IWP), INTENT(IN):: NDREF,NPREF, NG3
      Real(KIND=WP), INTENT(IN):: DREF(NDREF),PREF(NPREF),F3(NG3)
      Real(KIND=WP), INTENT(IN):: FD(NDREF),FP(NPREF)
      INTEGER(KIND=Byte), INTENT(IN):: idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      Real(KIND=WP) Dummy(1)
      INTEGER(KIND=IWP) MYRANK,MA
#endif
      INTEGER(KIND=IWP) ILO,IHI,JLO,JHI,LDA
      INTEGER(KIND=IWP) ICASE,ISYM,NIN,NAS,NBA,lg_BA,MBA
      Real(KIND=WP) DBA
      Real(KIND=WP), EXTERNAL:: PSBMAT_FPRINT

      ICASE=1
! LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NBA=(NAS*(NAS+1))/2
        IF(NBA.LE.0) CYCLE

! Set up the matrix BA(tuv,xyz) defined by the expression
! <ituv|H0-E0|kxyz> = dik ( alpha(a) SA(tuv,xyz) + BA(tuv,xyz) )
! Formula used:
! BA(tuv,xyz) = (Ey+Eu+Ex+Et-EASUM)*SA(tuv,xyz) - Fvuxtyz
! - dyu ( Fvzxt - Eu Gvzxt ) - dyt ( Fvuxz - Et Gvuxz )
! - dxu ( Fvtyz - Eu Gvtyz ) - dxu dyt ( Fvz - (Et+Eu) Gvz )
! + 2dxt ( Fvuyz - Et Gvuyz ) + 2dxt dyu ( Fvz - (Et+Eu) Gvz )

! where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
! Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

        CALL PSBMAT_GETMEM('BA',lg_BA,NAS)
        CALL PSBMAT_READ('S',iCase,iSym,lg_BA,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_BA,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(6,*) 'MKBA: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_BA,ILO,IHI,JLO,JHI,MA,LDA)
            MBA=LDA*(JHI-JLO+1)
            CALL MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,              &
     &                   DBL_MB(MA),MBA,                                &
     &                   ILO,IHI,JLO,JHI,LDA)
            CALL MKBA_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,NAS,LDA,           &
     &                       NG3,F3,IDXG3)
            CALL GA_RELEASE_UPDATE (LG_BA,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKBA_F3_MPP(ISYM,DUMMY,ILO,IHI,NAS,LDA,                &
     &                       NG3,F3,IDXG3)
          END IF
        ELSE
#endif
          LDA=0
          ILO=1
          IHI=NAS
          JLO=1
          JHI=NAS
          MBA=NAS*(NAS+1)/2
          CALL MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,                     &
     &                 ISYM,GA_Arrays(lg_BA)%A(:),MBA,                  &
     &                 ILO,IHI,JLO,JHI,LDA)
          CALL MKBA_F3(ISYM,GA_Arrays(lg_BA)%A(:),MBA,NG3,F3,IDXG3)
#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('B',iCase,iSYM,lg_BA,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DBA=PSBMAT_FPRINT(lg_BA,NAS)
          WRITE(6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'A', ISYM, DBA
        END IF

        CALL PSBMAT_FREEMEM(lg_BA)
      END DO

      END SUBROUTINE MKBA
