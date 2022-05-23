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
      SUBROUTINE OEISG(REL,SREL,TREL,UREL,ZETA,
     &                 ZN, MX100, NSYM, NBAS1,
     &                 unrel, tnrel, hcorr, iprint,
     &                  VEXTT,PVPT,EVN1,EVN2,RE1R,AUXI,W1W1,W1E0W1)

      IMPLICIT REAL*8 (A-H,O-Z)
C
      dimension REL(nBAS1*(nBAS1+1)/2),SREL(nBAS1*(nBAS1+1)/2),
     &          TREL(nBAS1*(nBAS1+1)/2),UREL(nBAS1*(nBAS1+1)/2),
     &          ZETA(nBAS1),
     &          unrel(nBAS1*(nBAS1+1)/2), tnrel(nBAS1*(nBAS1+1)/2),
     &          hcorr(nBAS1*(nBAS1+1)/2)
C
      dimension facto(16)
#include "relmp.fh"
      Real*8 VEXTT(*),PVPT(*),
     &       EVN1(NBAS1,NBAS1),EVN2(NBAS1,NBAS1),
     &       RE1R(NBAS1,NBAS1),
     &       AUXI(NBAS1,NBAS1),
     &       W1W1(NBAS1,NBAS1),W1E0W1(NBAS1,NBAS1)
C
C     ONE CONFIGURATION
C     ONE ELECTRON INTEGRAL PROGRAM
C     COMPUTES THE MATRICES S,U,T WITH SLATER (NFLAG1=0) OR GAUSSIAN
C      (NFLAG1=1) ORBITALS AS  BASIS
C
C
C
      CALL OEISG_INTERNAL(REL)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE OEISG_INTERNAL(REL)
      USE ISO_C_BINDING
      REAL*8, TARGET :: REL(*)
      INTEGER, POINTER :: iREL(:)
c     initialize
      CALL RELOP
c     coefficients
      FACTO(1)=1.D0
      FACTO(2)=1.D0
      DO 1301 I=3,16
      IM1=I-1
      FACTO(I)=IM1*FACTO(I-2)
 1301 CONTINUE
      iparm=2
*     WRITE (6,*) ' MX100=',MX100
*     write(6,*) ' nsym',nsym
      if(iprint.ge.10) then
         write(6,*) ' symmetry', nsym
         write(6,*) ' number of basis functions', nbas1
         write(6,*) ' charge', zn
         do i=1,nbas1
            write(6,*) zeta(i)
         enddo
      endif
      L=NSYM
C
      IJ=0
      np=L
      NQ=NP
      DO 10 I=1,NBAS1
      ZP=ZETA(I)
      WP=2.D0*(NP-L)/ZP
      DO 11 J=1,I
      IJ=IJ+1
      ZQ=ZETA(J)
      ZPQ=0.5D0*(ZP+ZQ)
      NPQ=NP+NQ+1
      WQ=2.D0*(NQ-L)/ZQ
      NPQ1=NPQ-1
      NPQ2=NPQ-2
      VPQ=FACTO(NPQ1)/ZPQ**(0.5D0*NPQ)
      VP=FACTO(2*NP)/ZP**(NP+0.5D0)
      VQ=FACTO(2*NQ)/ZQ**(NQ+0.5D0)
      IF (NPQ1.GT.2) GO TO 3
      VPQM2=1.D0
      GO TO 4
    3 VPQM2=FACTO(NPQ2-1)/ZPQ**(0.5D0*NPQ2)
    4 VPQ1=FACTO(NPQ2)*1.595769121605731D0/ZPQ**(0.5D0*NPQ1)
      VPQP2=FACTO(NPQ+1)/ZPQ**(0.5D0*(NPQ+2))
      TERM2=VPQP2-VPQ*(WP+WQ)+WP*WQ*VPQM2
      TERM1=1.D0/SQRT(VP*VQ)
      SREL(IJ)=TERM1*VPQ
      UREL(IJ)=TERM1*VPQ1
      UNREL(IJ)=TERM1*VPQ1
      TNREL(IJ)=0.5D0*ZP*ZQ*TERM1*TERM2
C.    TNRE=0.5D0*ZP*ZQ*TERM1*TERM2
      IF (IPARM.GT.0)   THEN
         I1=np-1
         I2=0
         I3=0
         J1=nq-1
         J2=0
         J3=0
C*
C*     PVP FOR XY-FUNCTIONS
C*
         if (np.eq.3) then
            REL(IJ)=EXTC(L,ZP,ZQ,1,1,0,1,1,0)
            GOTO 900
         endif
C*
C*     PVP FOR XYZ-FUNCTIONS
C*
         if (np.eq.4) then
            REL(IJ)=EXTC(L,ZP,ZQ,1,1,1,1,1,1)
            GOTO 900
         endif
C*
C*     PVP FOR X**(NP-1)-FUNCTIONS
C*
         REL(IJ)=EXTC(L,ZP,ZQ,I1,I2,I3,J1,J2,J3)
C.       TREL(IJ)=SQROPY(ZP,ZQ,I1,I2,I3,J1,J2,J3)
900      TREL(IJ)=0.5D0*ZP*ZQ*TERM1*TERM2
      ELSE
         TREL(IJ)=0.5D0*ZP*ZQ*TERM1*TERM2
C.       TNRE=0.5D0*ZP*ZQ*TERM1*TERM2
      ENDIF
   11 CONTINUE
   10 CONTINUE
C
C     CONVERT TO RELATIVISTIC INTEGRALS
C
      NBSZ=NBAS1*(NBAS1+1)/2
      NBSZ5=5*NBSZ
      NBSQ=NBAS1*NBAS1
      NBSQ5=5*NBSQ
      CALL C_F_POINTER(C_LOC(REL(NBSZ+1)),iREL,[NBSZ])
      CALL AT34R(NBAS1,NBSZ  ,ZN,
     *SREL,UREL,TREL,REL(1),
C     OVERLAP    POTENTIAL  KIN.ENERGY PVP
     *iREL,REL(2*NBSZ+1),REL(3*NBSZ+1),REL(4*NBSZ+1),
C     MULT      BU          P           G
     *REL(NBSZ5+1),REL(NBSZ5+NBSQ+1),REL(NBSZ5+2*NBSQ+1),
C     EIG            SINV                REVT
     *REL(NBSZ5+3*NBSQ+1),REL(NBSZ5+4*NBSQ+1),REL(NBSZ5+5*NBSQ+1),
C     AUX                   OVE                   EW
     *REL(NBSZ5+NBSQ5+MX100+1),REL(NBSZ5+NBSQ5+2*MX100+1),
C     E                         AA
     *REL(NBSZ5+NBSQ5+3*MX100+1),REL(NBSZ5+NBSQ5+4*MX100+1),iprint,
C    RR                         TT
     *VEXTT,PVPT,EVN1,EVN2,RE1R,AUXI,W1W1,W1E0W1)
      NULLIFY(iREL)
C
C     TRANSFER RELATIVISTIC KINETIC ENERGY AND POTENTIAL INTEGRALS
C
      if(iprint.ge.10) then
         write(6,*) ' matrices'
          write(6,*)  l,nbas1
         write(6,*) ' srel'
          write(6,12) (srel(j),j=1,(nbas1*(nbas1+1))/2)
         write(6,*) ' trel'
          write(6,12) (trel(j),j=1,(nbas1*(nbas1+1))/2)
         write(6,*) ' urel'
          write(6,12) (urel(j),j=1,(nbas1*(nbas1+1))/2)
         write(6,*) ' tnrel'
          write(6,12) (tnrel(j),j=1,(nbas1*(nbas1+1))/2)
         write(6,*) ' unrel'
          write(6,12) (unrel(j),j=1,(nbas1*(nbas1+1))/2)
         write(6,*) ' rel'
          write(6,12) (rel(j),j=1,(nbas1*(nbas1+1))/2)
      endif
      ij=0
      do 26 i=1,nbas1
      do 27 j=1,i
      ij=ij+1
      hcorr(ij)=(trel(ij)-zn*urel(ij)) - (tnrel(ij) - zn*unrel(ij))
      trel(ij)=trel(ij) - tnrel(ij)
      urel(ij)=-(zn*urel(ij) - zn*unrel(ij))
 27   continue
 26   continue
      if(iprint.ge.20) then
          write(6,*) ' full correction metrix'
          write(6,*)  l,nbas1
          write(6,12) (hcorr(j),j=1,(nbas1*(nbas1+1))/2)
      endif
 12   format(4f18.14)
      RETURN
      END SUBROUTINE OEISG_INTERNAL
*
      END
