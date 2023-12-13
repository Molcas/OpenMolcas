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
* Copyright (C) 2015, Kamal Sharkas                                    *
*               2019, Thomas J. Duignan                                *
*               2021, Rulin Feng                                       *
************************************************************************
*
* Note: The hyperfine code is based on the analogous
* pre-existing G-tensor functionality
*
      SUBROUTINE HFCTS(PROP,USOR,USOI,ENSOR,NSS,ENERGY,JBNUM,DIPSOM,
     &                 ESO,XYZCHR,BOLTZ_K)
      use rassi_aux, only: ipglob
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
      parameter (THRSH=1.0D-10)
      parameter (ZERO=0.0D0)
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "constants.fh"
#include "hfc_logical.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP),ENERGY(NSTATE),
     &          JBNUM(NSTATE)
      Character*1 xyzchr(3)
      Character*8 SDPROP
      Character*8 PSOPROP
      Character*8 DMPPROP
      Dimension IZMR(3),IZMI(3)
      Dimension GTENS(3,3)
      Dimension TMPMAT(3,3),TMPVEC(3,3),EVR(3),EVI(3)
      COMPLEX*16 ZEKL(2,2,3,NSTATE),GCONT(9,NSTATE)
      COMPLEX*16 DIPSOm(3,NSS,NSS),Z(NSS,NSS)
      COMPLEX*16 SPNSFS(3,NSS,NSS)
      COMPLEX*16 DIPSOf(3,NSS,NSS),DIMSO(3,3,NSS,NSS)
      COMPLEX*16 DIPSOfc(3,NSS,NSS),DIPSOfsd(3,NSS,NSS)
      COMPLEX*16 DIPSOfcsd(3,NSS,NSS),DIPSOfpso(3,NSS,NSS)
      REAL*8 GTOTAL(9),ESO(NSS)
      Dimension TMPf(NTP)
      Dimension HFC_1(3,3),HFC_2(3,3),HFC_3(3,3)
      Dimension CurieT(3,3),DiamT(3,3),PNMRCPS(NTP,NSS,3,3)
      Dimension PNMRT(NTP,3,3),PNMR(NTP,3,3)
      Dimension PNMRC(NTP,3,3),PNMRD(NTP,3,3)
      REAL*8 DLTTA,Zstat,p_Boltz,Boltz_k
      LOGICAL ISGS(NSS)
!      Dimension IMR(3),IMI(3)
      INTEGER IFUNCT
      REAL*8, Allocatable:: SOPRR(:,:), SOPRI(:,:)

!      AVOGADRO=CONST_AVOGADRO_
!      AU2EV=CONV_AU_TO_EV_
      AU2CM=CONV_AU_TO_CM1_
!      AU2T=CONV_AU_TO_T_
!      AU2J=CONV_AU_TO_KJ_*1.0D3
!      J2CM=AU2CM/AU2J
!      AU2JTM=(AU2J/AU2T)*AVOGADRO
      ALPHA=CONST_AU_VELOCITY_IN_SI_/CONST_C_IN_SI_
      ALPHA2= ALPHA*ALPHA
!      DEBYE=CONV_AU_TO_DEBYE_
!      AU2REDR=2.0D2*DEBYE
!      HALF=0.5D0

!      coeff_chi=0.1D0*AVOGADRO/CONST_BOLTZMANN_*
!     &          CONST_BOHR_MAGNETON_IN_SI_**2
      FEGVAL=-(CONST_ELECTRON_G_FACTOR_)
!      BOLTZ=CONST_BOLTZMANN_/AU2J
!      Rmu0=4.0D-7*CONST_PI_

      IF(IFSONCINI) THEN
      WRITE(6,*)
      WRITE(6,*) '  Soncini pNMR Tensor and A-Matrix Approach II '
      WRITE(6,*) '  ============================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ELSE
      WRITE(6,*)
      WRITE(6,*) '  A-Matrix Approach II                     '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ENDIF

C Mapping from spin states to spin-free state and to spin:
      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)
      ISS=0
      DO ISTATE=1,NSTATE
       JOB=JBNUM(ISTATE)
       MPLET=MLTPLT(JOB)
       DO MSPROJ=-MPLET+1,MPLET-1,2
        ISS=ISS+1
        IWORK(LMAPST-1+ISS)=ISTATE
        IWORK(LMAPSP-1+ISS)=MPLET
        IWORK(LMAPMS-1+ISS)=MSPROJ
       END DO
      END DO


      IAMX=0

      DO IPROP=1,NPROP
        IF(PNAME(IPROP)(1:3).EQ.'ASD'.AND.ICOMP(IPROP).EQ.1) THEN

* Get the center number
         Read (PNAME(IPROP),'(4x,i4)') ICEN

      WRITE(6,*) '  ========================================='
      write(6,*) '  A(Total)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='



C Identify which properties are ASD matrix elements:
c Labeled AMFI for now
c 1,2,3,4,5,6 -> xx,xy,xz,yy,yz,zz
      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
C Identify which properties are Orbital Paramagnetic (PSOP) matrix elements:
      IAMX=0
      IAMY=0
      IAMZ=0

      WRITE(SDPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",SDPROP
      WRITE(PSOPROP,'(a4,i4)') 'PSOP',ICEN
      WRITE(6,*) "Looking for ",PSOPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP)(1:3).EQ.SDPROP(1:3)
     &   .AND.PNAME(KPROP)(5:8).EQ.SDPROP(5:8)) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       ELSE IF(PNAME(KPROP).EQ.PSOPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMX=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMY=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMZ=KPROP
       END IF
      END DO

cccccccccccccccccccccccccccccccccccccccc
c Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc
c      SDPROP = 'MLTPL  0'
c      WRITE(6,*) "Looking for overlap matrix ",SDPROP
c      DO KPROP=1,NPROP
c       IF(PNAME(KPROP).EQ.SDPROP) THEN
c         IAMFI1=KPROP
c         IAMFI2=0
c         IAMFI3=0
c         IAMFI4=0
c         IAMFI5=0
c         IAMFI6=0
c       END IF
c      END DO
cccccccccccccccccccccccccccccccccccccccc
c End testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc


* These will hold the entire expression
* For g: L+2S
* For A: ?

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLZI),1)

      IF(IAMX.GT.0) CALL SMMAT(PROP,WORK(LLXI),NSS,IAMX,1)
      IF(IAMY.GT.0) CALL SMMAT(PROP,WORK(LLYI),NSS,IAMY,2)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,WORK(LLZI),NSS,IAMZ,3)


      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZI),1)


      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)


* Note that the 6th element contains r*r
* This is 4*pi * the contact term we want
* We need 8*pi/3 * < delta > so we just want
* a factor of 2/3
          ACNT = 0d0

          IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
           !write(6,*)"ACNT",ACNT
          ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
             WRITE(6,*) "********************"
             WRITE(6,*) "* Skipping FC Part *"
             WRITE(6,*) "********************"
          END IF

* Since this is a traceless tensor, generate the 6th
* element (3z^2-r^2)
          AMFI6=-AMFI1-AMFI4

          IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          ELSE
             IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
                WRITE(6,*) "********************"
                WRITE(6,*) "* Skipping SD Part *"
                WRITE(6,*) "********************"
             END IF
             AMFI1=0.0d0
             AMFI2=0.0d0
             AMFI3=0.0d0
             AMFI4=0.0d0
             AMFI5=0.0d0
             AMFI6=0.0d0
          END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          If (.not.MAG_X2C) Then
            CALL ADD_INFO("ASDFC1",[AMFI1],1,5)
            CALL ADD_INFO("ASDFC2",[AMFI2],1,5)
            CALL ADD_INFO("ASDFC3",[AMFI3],1,5)
            CALL ADD_INFO("ASDFC4",[AMFI4],1,5)
            CALL ADD_INFO("ASDFC5",[AMFI5],1,5)
            CALL ADD_INFO("ASDFC6",[AMFI6],1,5)
          End If

cccccccccccccccccccccccccccccccccccccccc
c Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc
c          AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
c          AMFI2=0.0d0
c          AMFI3=0.0d0
c          AMFI4=AMFI1
c          AMFI5=0.0d0
c          AMFI6=AMFI1
cccccccccccccccccccccccccccccccccccccccc
c END Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO


      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LZXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LZYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LZZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)


*     SVC 20090926 Experimental
*     Add analysis of different contributions

*     Establish which spin components of SFS belong to the ground state
      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

*     Analyze the different contributions to the GS Kramers doublet
*     Zeeman matrix elements.  There are 4 different ME's: <1|Ze|1>,
*     <1|Ze|2>, <2|Ze|1>, and <2|Ze|2>, stored in ZEKL.  Contributions
*     of SFS i,j to SOS k,l (k,l=1,2): <k|Ze|l> = Sum(i,j) U(i,k)*
*     <i|Ze|j> U(j,l).  This sum is decomposed into parts belonging to
*     each SFS state i as follows:
*     -> GS's contain only MEs with themselves and other GS's
*     -> ES's contain MEs with themselves, the GS's (2x) and other ES's
*        The ME's with the GS's are counted twice as they do not belong to
*        any GS's (they contain only ME's within their own GS group)
*        The contributions with other ES's are split between the ES's,
*        counting them double (<i|Ze|j> and <j|Ze|i>) and divide by two later.

      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C       DO IXYZ=1,3
C        WRITE(6,FMT=720) 'ZEKL', IXYZ, ISTATE,
C     &      ZEKL(:,:,IXYZ,ISTATE)
C       ENDDO

C 710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
C 720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

*     We now have decomposed the <k|Ze|l> into terms belonging to either
*     a GS or an ES for each k,l=1,2 and p=x,y,z stored in ZEKL(k,l,p,SFS)
*     Now, these new decomposed terms of the ZEKL ME's are combined to
*     form the G tensor.  Consider e.g. that <k|Ze|l> is decomposed into
*     <k|GS|l> + <k|ES1|l> + <k|ES2|l>, then the contributions to G are given as:
*     -> G_pq/2 = <k|Ze_p|l> <l|Ze_q|k>
*             = (<k|GS_p|l> + <k|ES1_p|l> + <k|ES2_p|l>)
*             * (<l|GS_q|k> + <l|ES1_q|k> + <l|ES2_q|k>)
*
*     from GS: (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2
*     from ES1: 2*((<k|ES1_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_q|l>/2 * <l|ES1_p|k>/2)/2)
*               + (<k|ES1_p|l>/2 * <l|ES1_q|k>/2)/2 + (<k|ES1_q|l>/2 * <l|ES2_p|k>/2)/2
*               + (<k|ES1_p|l>/2 * <l|ES2_q|k>/2)/2 + (<k|ES2_q|l>/2 * <l|ES1_p|k>/2)/2
*     In the end, the outer division by 2 cancels on both sides, and the
*     inner divisions by two combine to a division by 4.

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE

* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOf(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo

c
c*     Continue original calculation of G tensor (=gg^*)
c

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL PRCMAT(NSS,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      IF(IFSONCINI) THEN

      do ic=1,3
      do jc=1,3
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIMSO(ic,jc,ISS,JSS)=0.0d0
      enddo
      enddo
      enddo
      enddo
      WRITE(DMPPROP,'(a4,i4)') 'DMP ',ICEN
      WRITE(6,*) "Looking for ",DMPPROP
      DO KPROP=1,NPROP
      IF(PNAME(KPROP).EQ.DMPPROP) THEN
        Call mma_Allocate(SOPRR,NSS**2,1,Label='SOPRR')
        Call mma_Allocate(SOPRI,NSS**2,1,Label='SOPRI')
        SOPRR(:,1)=0.0D0
        SOPRI(:,1)=0.0D0

      CALL SMMAT(PROP,SOPRR,NSS,KPROP,0)
      CALL ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
      !CALL PRCMAT(NSS,SOPRR,SOPRI)
      CALL MULMAT(NSS,SOPRR,SOPRI,eex,Z)
      ic=0
      jc=0
      !write(6,*)"ICOMP(KPROP)",ICOMP(KPROP)
      if (ICOMP(KPROP).EQ.1) then
      ic=1
      jc=1
      elseif (ICOMP(KPROP).EQ.2) then
      ic=1
      jc=2
      elseif (ICOMP(KPROP).EQ.3) then
      ic=1
      jc=3
      elseif (ICOMP(KPROP).EQ.4) then
      ic=2
      jc=1
      elseif (ICOMP(KPROP).EQ.5) then
      ic=2
      jc=2
      elseif (ICOMP(KPROP).EQ.6) then
      ic=2
      jc=3
      elseif (ICOMP(KPROP).EQ.7) then
      ic=3
      jc=1
      elseif (ICOMP(KPROP).EQ.8) then
      ic=3
      jc=2
      else
      ic=3
      jc=3
      endif
      DO ISS=1,NSS
      DO JSS=1,NSS

      if (ISS.eq.JSS) then
      !DIMSO(ICOMP(KPROP),ISS,JSS)=Z(ISS,JSS)
      DIMSO(ic,jc,ISS,JSS)=Z(ISS,JSS)
      ELSE
      DIMSO(ic,jc,ISS,JSS)=0
      ENDIF
      enddo
      enddo
      CAll mma_deallocate(SOPRR)
      CAll mma_deallocate(SOPRI)
      ENDIF

      END DO
      iT=0
      do iT=1,NTP
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=0.d0
      PNMR(iT,ic,jc)=0.d0
      PNMRC(iT,ic,jc)=0.d0
      PNMRD(iT,ic,jc)=0.d0
      enddo
      enddo
      enddo

      do iT=1,NTP
      do Iss=1,Nss
      do ic=1,3
      do jc=1,3
      PNMRCPS(iT,Iss,ic,jc)=0.d0
      enddo
      enddo
      enddo
      enddo

      iT=0
      do iT=1,NTP
      if(iT.eq.1) then
      TMPf(iT)=TMINP+0.0001d0
      ELSE
      DLTTA=(TMAXP-TMINP)/(dble(NTP-1))
      TMPf(iT)=TMINP+DLTTA*dble(iT-1)
      ENDIF
      Zstat=0.d0
      do Iss=1,Nss
      p_Boltz=EXP(-ESO(Iss)/Boltz_k/TMPf(iT))
      !write(6,*)'p_Boltz',p_Boltz
      Zstat=Zstat+p_Boltz
      do IC=1,3
      do JC=1,3
      DiamT(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      CurieT(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      HFC_2(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      HFC_3(IC,JC)=0.d0
      enddo
      enddo
      do Jss=1,Nss
      dlt_E=Eso(Iss)-Eso(Jss)
      do IC=1,3
      do JC=1,3
      HFC_1(IC,JC)=0.d0
      enddo
      enddo
      do ic=1,3
      do jc=1,3

      HFC_1(ic,jc)=DBLE(DIPSOm(ic,Iss,Jss)*
     & CONJG(DIPSOf(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))

!!      HFC_1(ic,jc)=DBLE(DIPSOf(ic,Iss,Jss)*
!!     & CONJG(DIPSOm(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))


!      HFC_3(ic,jc)= 1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))/(ALPHA2*AU2CM)

      !if(ABS(dlt_E).LT.1.D-3) then
      if(ABS(dlt_E).LT.10.97D0) then
      CurieT(ic,jc)=   CurieT(ic,jc)+HFC_1(ic,jc)
      HFC_2(ic,jc)=    HFC_2(ic,jc)+ HFC_1(ic,jc)
      HFC_3(ic,jc)=    HFC_3(ic,jc)+ 0.d0*HFC_1(ic,jc)
      DiamT(ic,jc)=    DiamT(ic,jc)+ 0.d0*HFC_1(ic,jc)
      else

      HFC_2(ic,jc)= HFC_2(ic,jc)-(DBLE(CONJG(DIPSOm(ic,Iss,Jss))*
     & DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & CONJG(DIPSOf(jc,Iss,Jss))))/dlt_E
!      HFC_2(ic,jc)= HFC_2(ic,jc)-2.d0*(DBLE(CONJG(DIPSOm(ic,Iss,Jss))*
!     & DIPSOf(jc,Iss,Jss)))/dlt_E

      HFC_3(ic,jc)= HFC_3(ic,jc)-(DBLE(CONJG(DIPSOm(ic,Iss,Jss))*
     & DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & CONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      CurieT(ic,jc)= CurieT(ic,jc)-0.d0*
     &(DBLE(CONJG(DIPSOm(ic,Iss,Jss))*
     &DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & CONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      DiamT(ic,jc)=DiamT(ic,jc)-0.d0*
     &(DBLE(CONJG(DIPSOm(ic,Iss,Jss))*
     &DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & CONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      endif

      HFC_2(ic,jc)= HFC_2(ic,jc) + (1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))
     & /(ALPHA2*AU2CM))
      DiamT(ic,jc)=DiamT(ic,jc)+ (1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))
     & /(ALPHA2*AU2CM))
      enddo
      enddo
      enddo !Jss
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=    PNMRT(iT,ic,jc)+ p_Boltz*
     & HFC_2(ic,jc)

      PNMR(iT,ic,jc)=    PNMR(iT,ic,jc)+ p_Boltz*
     & HFC_3(ic,jc)

      PNMRC(iT,ic,jc)=    PNMRC(iT,ic,jc)+ p_Boltz*
     & CurieT(ic,jc)

      PNMRCPS(iT,Iss,ic,jc)= PNMRCPS(iT,Iss,ic,jc)+CurieT(ic,jc)
      PNMRCPS(iT,Iss,ic,jc)=1.D6*AU2CM*ALPHA2*PNMRCPS(iT,Iss,ic,jc)

      PNMRD(iT,ic,jc)=    PNMRD(iT,ic,jc)+ p_Boltz*
     & DiamT(ic,jc)

      enddo
      enddo
      enddo !Iss
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRT(iT,ic,jc)/Zstat)
      PNMR(iT,ic,jc) =1.D6*AU2CM*ALPHA2*(PNMR(iT,ic,jc)/Zstat)
      PNMRC(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRC(iT,ic,jc)/Zstat)
      PNMRD(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRD(iT,ic,jc)/Zstat)
      !write(6,*) PNMRT(iT,ic,jc)
      !AU2CM*ALPHA2*
      enddo
      enddo

      enddo ! iT

      !IT=0
      !IC=0
      !JC=0
      !Do IT=1,NTP
      !Do IC=1,3
      !Do JC=1,3
      !PNMRT(IT,IC,JC)=1.d6*PNMRT(IT,IC,JC)
      !enddo
      !enddo
      !enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') ' TOTAL SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),
     & ((PNMRT(iT,ic,jc),jc=1,3),ic=1,3)
      enddo
      !write(6,*)'*******'
      !do iT=1,NTP
      !do ic=1,3
      !do jc=1,3
      !write(6,*)PNMRT_tens(iT,ic,jc)
      !enddo
      !enddo
      !enddo
      !write(6,*)'*******'

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'paramagnetic SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),
     & ((PNMR(iT,ic,jc),jc=1,3),ic=1,3)
      enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Diamagnetic  SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),
     & ((PNMRD(iT,ic,jc),jc=1,3),ic=1,3)
      enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Curie term contrib. to SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,9X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),
     & ((PNMRC(iT,ic,jc),jc=1,3),ic=1,3)
!       IF (.EQV..TRUE.) THEN
      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Curie per state in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'NSS','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      !write(6,*)""
      !write(6,*) PNMRC(iT,3,3)
      !write(6,*)""
      do Iss=1,NSS
      write(6,'(4X,I5,3X,11(F10.4,2X),F8.4)') Iss,
     & ((PNMRCPS(iT,Iss,ic,jc),jc=1,3),ic=1,3)
       enddo

!       ENDIF
       write(6,'(/)')
       write(6,'(10A)') (('============='), K=1,10)
       write(6,'(/)')
      enddo
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '  A-Matrix'
      WRITE(6,*) '  =========='
      ENDIF !IFSONCINI

! 1111 Continue

      ISS=1
      DO WHILE ( (ISS.LE.NSS) .AND.
     &          ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRATHR) )

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      !WRITE(6,*) 'KDGN=',KDGN
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFATCALSA) GOTO 451
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOf,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  451 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 1780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF


      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

C construct g_s matrix from G by back-transormation of the
C square root of the G eigenvalues
      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

      CALL ADD_INFO("ATENS",GTENS,9,5)
      CALL ADD_INFO("ATENS2",EVR,3,5)

 1780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

      IF(.NOT.IFACALFCON) GOTO 1901

      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (FC)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(SDPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",SDPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP)(1:3).EQ.SDPROP(1:3)
     &   .AND.PNAME(KPROP)(5:8).EQ.SDPROP(5:8)) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       END IF
      END DO

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLZI),1)

      IF(IAMX.GT.0) CALL SMMAT(PROP,WORK(LLXI),NSS,IAMX,0)
      IF(IAMY.GT.0) CALL SMMAT(PROP,WORK(LLYI),NSS,IAMY,0)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,WORK(LLZI),NSS,IAMZ,0)

      CALL GETMEM('MXR','ALLO','REAL',LMXR,NSS**2)
      CALL GETMEM('MXI','ALLO','REAL',LMXI,NSS**2)
      CALL GETMEM('MYR','ALLO','REAL',LMYR,NSS**2)
      CALL GETMEM('MYI','ALLO','REAL',LMYI,NSS**2)
      CALL GETMEM('MZR','ALLO','REAL',LMZR,NSS**2)
      CALL GETMEM('MZI','ALLO','REAL',LMZI,NSS**2)

!      IMR(1)=LMXR
!      IMI(1)=LMXI
!      IMR(2)=LMYR
!      IMI(2)=LMYI
!      IMR(3)=LMZR
!      IMI(3)=LMZI

      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMXI),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMYI),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LMZI),1)

      CALL SMMAT(PROP,WORK(LMXR),NSS,0,1)
      CALL SMMAT(PROP,WORK(LMYI),NSS,0,2)
      CALL SMMAT(PROP,WORK(LMZR),NSS,0,3)

      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMXR),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMYI),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMZR),1)

      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LMXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LMYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LMZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMXR),WORK(LMXI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMYR),WORK(LMYI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMZR),WORK(LMZI))

      CALL GETMEM('MXR','FREE','REAL',LMXR,NSS**2)
      CALL GETMEM('MXI','FREE','REAL',LMXI,NSS**2)
      CALL GETMEM('MYR','FREE','REAL',LMYR,NSS**2)
      CALL GETMEM('MYI','FREE','REAL',LMYI,NSS**2)
      CALL GETMEM('MZR','FREE','REAL',LMZR,NSS**2)
      CALL GETMEM('MZI','FREE','REAL',LMZI,NSS**2)

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0
* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
         ! END IF
          !AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
           !   AMFI1=-1.0*AMFI1
           !   AMFI2=-1.0*AMFI2
           !   AMFI3=-1.0*AMFI3
           !   AMFI4=-1.0*AMFI4
           !   AMFI5=-1.0*AMFI5
           !   AMFI6=-1.0*AMFI6
          !ELSE
          !   IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !      WRITE(6,*) "********************"
          !      WRITE(6,*) "* Skipping SD Part *"
          !      WRITE(6,*) "********************"
          !   END IF
             AMFI1=0.0d0
             AMFI2=0.0d0
             AMFI3=0.0d0
             AMFI4=0.0d0
             AMFI5=0.0d0
             AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO

      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO
      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF

C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C       DO IXYZ=1,3
C        WRITE(6,FMT=720) 'ZEKL', IXYZ, ISTATE,
C     &      ZEKL(:,:,IXYZ,ISTATE)
C       ENDDO

C 710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
C 720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE

* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfc(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRATHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      !WRITE(6,*) 'KDGN=',KDGN
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFATCALSA) GOTO 452
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfc,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  452 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 2780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF


      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
       GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &  TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO


 2780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1901  CONTINUE

      IF(.NOT.IFACALSDON) GOTO 1902

      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (SD)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(SDPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",SDPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP)(1:3).EQ.SDPROP(1:3)
     &   .AND.PNAME(KPROP)(5:8).EQ.SDPROP(5:8)) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       END IF
      END DO

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
          !   ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
          !END IF
          AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          !ELSE
           !  IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
            !    WRITE(6,*) "********************"
             !   WRITE(6,*) "* Skipping SD Part *"
              !  WRITE(6,*) "********************"
             !END IF
             !AMFI1=0.0d0
             !AMFI2=0.0d0
             !AMFI3=0.0d0
             !AMFI4=0.0d0
             !AMFI5=0.0d0
             !AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO

      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF

C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C       DO IXYZ=1,3
C        WRITE(6,FMT=720) 'ZEKL', IXYZ, ISTATE,
C     &      ZEKL(:,:,IXYZ,ISTATE)
C       ENDDO

C 710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
C 720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

        IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE

* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfsd(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRATHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFATCALSA) GOTO 453
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfsd,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  453 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 3780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF


      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
       GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      !DO IXYZ=1,3
      !   write(6,*) 'IXYZ',IXYZ
      !   write(6,*) ''
      ! DO JXYZ=1,3
      !   write(6,*) 'JXYZ',JXYZ
      !   write(6,*) ''
      ! GTENS(IXYZ,JXYZ)=0.0D0
      !  DO KXYZ=1,3
      !   write(6,*) 'KXYZ',KXYZ
      !   write(6,*) ''
      !   write(6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
      !   write(6,*) ''
      !   write(6,*) 'TMPVEC(IXYZ,KXYZ)',TMPVEC(IXYZ,KXYZ),IXYZ,KXYZ
      !   write(6,*) ''
      !   write(6,*) 'SQRT(EVR(KXYZ))',SQRT(EVR(KXYZ)),KXYZ
      !   write(6,*) ''
      !   write(6,*) 'TMPVEC(JXYZ,KXYZ)',TMPVEC(JXYZ,KXYZ),JXYZ,KXYZ
      !   write(6,*) '*************************'
      !  GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
      !&   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
      !  write(6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
      !   write(6,*) '*************************'
      !  END DO
      ! END DO
      !ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

 3780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1902  CONTINUE

* Skip if not a hyperfine calculation
      IF(.NOT.IFACALFCSDON) GOTO 1903


      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (FC+SD)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(SDPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",SDPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP)(1:3).EQ.SDPROP(1:3)
     &   .AND.PNAME(KPROP)(5:8).EQ.SDPROP(5:8)) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
      END IF
      END DO

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
          !END IF

          AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          !ELSE
           !  IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
           !     WRITE(6,*) "********************"
           !     WRITE(6,*) "* Skipping SD Part *"
           !     WRITE(6,*) "********************"
           !  END IF
           !  AMFI1=0.0d0
           !  AMFI2=0.0d0
           !  AMFI3=0.0d0
           !  AMFI4=0.0d0
           !  AMFI5=0.0d0
           !  AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)


          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5

        END DO
      END DO

      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO
      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE
* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO


        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfcsd(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo


      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRATHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      !WRITE(6,*) 'KDGN=',KDGN
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFATCALSA) GOTO 454
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfcsd,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  454 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 4780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF
      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO
      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

C construct g_s matrix from G by back-transormation of the
C square root of the G eigenvalues
      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

 4780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

 1903  CONTINUE

* Skip if not a hyperfine calculation
      IF(.NOT.IFACALPSO) GOTO 1904


      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (PSO)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

       IAMX=0
       IAMY=0
       IAMZ=0

      WRITE(PSOPROP,'(a4,i4)') 'PSOP',ICEN
      WRITE(6,*) "Looking for ",PSOPROP
      DO KPROP=1,NPROP
         IF(PNAME(KPROP).EQ.PSOPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMX=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMY=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMZ=KPROP
       END IF
      END DO

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLZI),1)

      IF(IAMX.GT.0) CALL SMMAT(PROP,WORK(LLXI),NSS,IAMX,1)
      IF(IAMY.GT.0) CALL SMMAT(PROP,WORK(LLYI),NSS,IAMY,2)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,WORK(LLZI),NSS,IAMZ,3)

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LZZI),1)


      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LZXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LZYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LZZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO
******************************************************
      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C       DO IXYZ=1,3
C        WRITE(6,FMT=720) 'ZEKL', IXYZ, ISTATE,
C     &      ZEKL(:,:,IXYZ,ISTATE)
C       ENDDO

C 710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
C 720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO


      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=CMPLX(0.0d0,0.0d0,kind=8)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE

* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfpso(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo


      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRATHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      !WRITE(6,*) 'KDGN=',KDGN
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFATCALSA) GOTO 455
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfpso,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  455 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 5780
      ENDIF


      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF


      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)


        GTIJ=GTIJ+CONTRIB


        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

!      IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
!      END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

C construct g_s matrix from G by back-transormation of the
C square root of the G eigenvalues

      WRITE(6,*)

      ICOUNT = 0
      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        If (EVR(KXYZ).LT.THRSH) THEN ! EVR(KXYZ) should probably be ZERO
          IF (EVR(KXYZ).LT.-THRSH) THEN
           write(6,*) "WARNNING : negative real eigenvalue of G"
           write(6,*) "WARNNING : this may cause numerical problem"
           write(6,*) "WARNNING : and may give 'NAN' in G tensor"
          ELSE
           EVR(KXYZ)=ZERO
           If (ICOUNT.EQ.0) then
           write(6,*) " The eigenvalue 'EVR(KXYZ)' will be set to ZERO"
           write(6,*) " each time it is very small"
           write(6,*) " to avoid the roots of negative eigenvalues"
           endif
           ICOUNT=1
          ENDIF
        ENDIF
        GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO
      WRITE(6,*)
      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

      CALL ADD_INFO('ATENS_PSO',GTENS,9,5)
      CALL ADD_INFO('ATENS_PSOEVR',EVR,3,5)

       WRITE(6,*)
 5780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1904  CONTINUE

* End loop over CNT properties
      END IF
      END DO

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      END
