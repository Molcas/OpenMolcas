!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE dampF(r,RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)
      IMPLICIT NONE
      INTEGER NCMM,MMLR(NCMM),IVSR,IDSTT,IDFF,FIRST,                    &
     &  m,MM
      REAL*8 r,RHOAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),br,XP,YP,            &
     &  DM(NCMM),                                                       &
     &  bpm(20,-2:0), cpm(20,-2:0),ZK
       DATA bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
!     WRITE(6,*) 'Made it inside of dampF! IVSR=',IVSR
      IF(NCMM.GT.4) THEN
       WRITE(6,*) 'IDSTT=',IDSTT
      ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
!                     bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      bpm(m,IDFF)= bDS(IDFF)/DBLE(m)
!                     cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DBLE(m))
                  ENDDO
              ENDDO
              FIRST= 0
          ENDIF
          br= RHOAB*r
!         WRITE(6,*) 'NCMM=',NCMM
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IVSR) + cpm(MM,IVSR)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
!... Actually ...  DM(m)= YP**(MM + IVSR/2)  :  set it up this way to
!   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IVSR.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IVSR.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IVSR.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IVSR.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(IVSR.EQ.-9) THEN
                  ENDIF
          ENDDO
          br=bTT(1) !Make sure that it's "referenced" in subroutine too!
      RETURN
! 600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IVSR=',i3,' no dampi
!    1ng function is defined')
! 602 FORMAT( /,' ***ERROR ***  RHOAB=', F7.4,'  yields an invalid Dampi
!    1ng Function definition')
      END SUBROUTINE dampF
