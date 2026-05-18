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
! Case C (ICASE=4)
!*******************************************************************************
      SUBROUTINE MKBC(DREF,NDREF,PREF,NPREF,FD,FP,NG3,F3,idxG3)
      use definitions, only: iwp, wp, Byte
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM,NINDEP,NTUV
      IMPLICIT NONE
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER(KIND=IWP), INTENT(IN):: NDREF, NPREF, NG3
      Real(KIND=WP), INTENT(IN):: DREF(NDREF),PREF(NPREF),F3(NG3)
      Real(KIND=WP), INTENT(IN):: FD(NDREF),FP(NPREF)
      INTEGER(KIND=Byte), INTENT(IN):: idxG3(6,NG3)
#ifdef _MOLCAS_MPP_
      Real(KIND=WP) Dummy(1)
      INTEGER(KIND=IWP) MYRANK,MA
#endif
      INTEGER(KIND=IWP) ILO,IHI,JLO,JHI,LDA,MBC
      INTEGER(KIND=IWP) ICASE,ISYM,NIN,NAS,NBC,lg_BC
      Real(KIND=WP) DBC
      Real(KIND=WP), EXTERNAL:: PSBMAT_FPRINT

      ICASE=4
! LONG loop over superindex symmetry.
      DO ISYM=1,NSYM
        NIN=NINDEP(ISYM,ICASE)
        IF(NIN.EQ.0) CYCLE
        NAS=NTUV(ISYM)
        NBC=(NAS*(NAS+1))/2
        IF(NBC.LE.0) CYCLE

! Set up the matrix BC(tuv,xyz) defined by the expression
! <atuv|H0-E0|cxyz> = dac ( alpha(a) SC(tuv,xyz) + BC(tuv,xyz) )
! Formula used:
!    BC(tuv,xyz)
!    = Fvutxyz +dyu Fvztx + dyx Fvutz + dtu Fvxyz + dtu dyx Fvz
!    +(Ey+Eu-EASUM)*SC(tuv,xyz)
!    -Eu*( dyu Gvztx + dtu Gvxyz )
!    -Ey dyx Gvutz
!    -(Eu+Ey)*( dtu dyx Gvz )

! where dyu = Kronecker(y,u) etc. Gvutxyz=<Evutxyz>, etc.
! Similarly, Fvutxyz= Sum(w)(EPSA(w)<Evutxyzww>, etc.

        CALL PSBMAT_GETMEM('BC',lg_BC,NAS)
        CALL PSBMAT_READ('S',iCase,iSym,lg_BC,NAS)

        ! fill in the 3-el parts
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          MYRANK = GA_NODEID()
          CALL GA_DISTRIBUTION (LG_BC,MYRANK,ILO,IHI,JLO,JHI)
          IF (JLO.NE.0 .AND. (JHI-JLO+1).NE.NAS) THEN
            WRITE(6,*) 'MKBC: MISMATCH IN RANGE OF THE SUPERINDICES'
            CALL ABEND()
          END IF
          IF (ILO.GT.0 .AND. JLO.GT.0) THEN
            CALL GA_ACCESS (LG_BC,ILO,IHI,JLO,JHI,MA,LDA)
            MBC=LDA*(JHI-JLO+1)
            CALL MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,              &
     &                   DBL_MB(MA),MBC,                                &
     &                   ILO,IHI,JLO,JHI,LDA)
            CALL MKBC_F3_MPP(ISYM,DBL_MB(MA),ILO,IHI,NAS,LDA,           &
     &                       NG3,F3,IDXG3)
            CALL GA_RELEASE_UPDATE (LG_BC,ILO,IHI,JLO,JHI)
          ELSE
            CALL MKBC_F3_MPP(ISYM,DUMMY,ILO,IHI,NAS,LDA,                &
     &                       NG3,F3,IDXG3)
          END IF
        ELSE
#endif
          ILO=1
          IHI=NAS
          JLO=1
          JHI=NAS
          LDA=0
          MBC=NAS*(NAS+1)/2
          CALL MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,                     &
     &                 ISYM,GA_Arrays(lg_BC)%A(:),MBC,                  &
     &                 ILO,IHI,JLO,JHI,LDA)
          CALL MKBC_F3(ISYM,GA_Arrays(lg_BC)%A(:),MBC,NG3,F3,IDXG3)

#ifdef _MOLCAS_MPP_
        END IF
#endif

        CALL PSBMAT_WRITE('B',iCase,iSYM,lg_BC,NAS)

        IF(IPRGLB.GE.DEBUG) THEN
          DBC=PSBMAT_FPRINT(lg_BC,NAS)
          WRITE(6,'("DEBUG> ",A4,1X,I3,1X,ES21.14)') 'C', ISYM, DBC
        END IF

        CALL PSBMAT_FREEMEM(lg_BC)
      END DO

      END SUBROUTINE MKBC
