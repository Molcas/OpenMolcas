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
      SUBROUTINE LEVXPC(KV,JR,EPR,GAMA,NPP,WF,RFN,V,VLIM,YH,DREF,       &
     &                             NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Calculates expectation values of the kinetic energy and of X**IP
!  (IP=1,MORDR), denoted XPTKE and XPCTR(IP), respectively, for level
!  v=KV, J=JR, E=EPR(cm-1), using wave function WF(i), (i=NBEG,NEND).
!** Assumes units of length are (Angstroms) .
!** Division by BFCT converts potential V(I) to units (cm-1).
!** If (|LXPCT| = 2  or  4), "punch" (WRITE(7,XXX)) results to channel-7
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
      REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),SDRDY(:),    &
     & VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
!!!!!
      INTEGER I,K,IRFN,IPNCH,ITRY,JR,KV,LXPCT,NPP,NBEG,NEND,MORDR
      REAL*8  WF(NPP),RFN(NPP),V(NPP),XPCTR(0:11),DM(0:20)
      REAL*8 BFCT,DS,DRT,DMR,DER,EPR,EINN,GAMA,YH,DREF,                 &
     &  RR,RXPCT,SS2,SF2,VLIM,XPTKE,PINV
!
      NDIMR= 131074
      CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
      CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
      CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
      CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
      CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
      CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
      EINN= BFCT*EPR
      IPNCH=0
      IF((IABS(LXPCT).EQ.2).OR.(IABS(LXPCT).GE.4)) IPNCH=1
!** MORDR is the highest-power expectation value considered.
      IF(MORDR.GT.11) MORDR=11
      ITRY=20
      IF(((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9)))                &
     &                                    .AND. (DREF.LE.0.D0)) ITRY=0
!** Start by calculating contributions at end points
    2 SS2=WF(NBEG)**2 *DRDY2(NBEG)
      SF2=WF(NEND)**2 *DRDY2(NEND)
      XPTKE= 0.5D0*(SS2*(EINN-V(NBEG)) + SF2*(EINN-V(NEND)))
      IF(MORDR.GT.0) THEN
          XPCTR(0)= 1.d0/YH
          DO  K=1,MORDR
              SS2=SS2*RFN(NBEG)
              SF2=SF2*RFN(NEND)
              XPCTR(K)=0.5D0*(SS2+SF2)
              ENDDO
          ENDIF
      IF(IRFN.GT.-4) THEN
!** For regular expectation values of a radial function ...
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2 *DRDY2(I)
              XPTKE= XPTKE+ DS*(EINN-V(I))
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
        ELSE
!** For expectation values involving partial derivative operator ...
          DO  K= 0,MORDR
              XPCTR(K)= 0.d0
              ENDDO
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2 *DRDY2(I)
              XPTKE= XPTKE+ DS*(EINN-V(I))
              DS= WF(I)*(WF(I+1)- WF(I-1))*DRDY2(I)
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
          DO  K= 0,MORDR
              XPCTR(K)= XPCTR(K)/(2.d0*YH)
              ENDDO
        ENDIF
      XPTKE= XPTKE*YH/BFCT
      IF(MORDR.LT.0) GO TO 99
      DMR= 0.d0
      DO  K=0,MORDR
          XPCTR(K)=XPCTR(K)*YH
          DMR= DMR+ DM(K)*XPCTR(K)
          ENDDO
      IF((LXPCT.EQ.1).OR.(IABS(LXPCT).EQ.2)) THEN
          IF(EPR.LE.VLIM) WRITE(6,600) KV,JR,EPR,DMR,XPTKE
          IF(EPR.GT.VLIM) WRITE(6,602) KV,JR,EPR,DMR,XPTKE,GAMA
          IF(IABS(IRFN).LE.9) WRITE(6,604) (K,XPCTR(K),K=1,MORDR)
          IF(IPNCH.GE.1) WRITE(7,701) KV,JR,EPR,GAMA,XPTKE,DMR,         &
     &                                        (XPCTR(K),K=1,MORDR)
          ENDIF
      IF(ITRY.GT.19) GO TO 99
!** If appropriate, iteratively correct DREF value till distance
!  coordinate expectation value is identically zero.
      IF(IRFN.EQ.-1) THEN
!** For Dunham expansion parameter, define revised function here
          DREF=XPCTR(1)
          DRT=DREF
          WRITE(6,603) ITRY,DRT,DREF
          DO  I= 1,NPP
              RVB(I)= RVB(I)/DREF - 1.D0
              ENDDO
          ITRY=99
          GO TO 2
          ENDIF
!** For Surkus-type expansion parameter, define revised function
      ITRY=ITRY+1
      IF(ITRY.EQ.1) THEN
          RXPCT= XPCTR(1)
          DREF= 0.D0
          DRT= RXPCT
        ELSE
          DER= -IRFN/(2.d0*DREF)
          DRT= -XPCTR(1)/DER
        ENDIF
      DREF=DREF+DRT
      WRITE(6,603) ITRY,DRT,DREF
!** Redefine Surkus-type distance variable RFN using new DREF
      PINV= 1.d0/PRV
      WRITE(6,*) PINV ! Make sure it's "referneced" in THIS subroutine.
      DO  I= 1,NPP
          RFN(I)= (RVB(I)**IRFN - DREF**IRFN)/(RVB(I)**IRFN+ DREF**IRFN)
          ENDDO
      IF(DABS(DRT/DREF).GE.1.D-12) GO TO 2
      CALL MMA_DEALLOCATE(RVB)
      CALL MMA_DEALLOCATE(YVB)
      CALL MMA_DEALLOCATE(DRDY2)
      CALL MMA_DEALLOCATE(FAS)
      CALL MMA_DEALLOCATE(SDRDY)
      CALL MMA_DEALLOCATE(VBZ)
   99 RETURN
  600 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,       &
     &  '   <KE>=',F11.3)
  602 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,       &
     & '   <KE>=',F11.3/'   Tunneling predissociation  Width(FWHM)=',   &
     & G13.6,'    <X**',I2,'>=',F13.8)
  604 FORMAT((8x,3('   <X**',I2,'>=',F13.8:)))
  603 FORMAT(' On iteration #',I2,'  change DREF by',1PD10.2,           &
     &  '  to   DREF=',0PF13.10,' [Angstroms]')
  701 FORMAT(2I4,F11.3,G11.4,F11.3,3(F12.7)/(5X,6F12.7))
      END
