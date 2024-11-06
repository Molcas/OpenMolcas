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
      SUBROUTINE SONATORBM_INT(DENS, CHARPROP, IC_,CHARTYPE,ASS,BSS,
     &                         iOpt,ROTMAT,
     &                         PROPVALXR,PROPVALYR,PROPVALZR,
     &                         PROPVALXI,PROPVALYI,PROPVALZI)
      use OneDat, only: sOpSiz
      use stdalloc, only: mma_allocate, mma_deallocate

      IMPLICIT None
#include "rassi.fh"
      Real*8 DENS(6,NBTRI)
      CHARACTER(LEN=8) CHARPROP
      INTEGER IC_
      CHARACTER(LEN=8) CHARTYPE
      INTEGER ASS, BSS, iOpt
      Real*8 ROTMAT(3,3),
     &       PROPVALXR,PROPVALYR,PROPVALZR,
     &       PROPVALXI,PROPVALYI,PROPVALZI

      Integer IDUM(1)
      Real*8, allocatable:: IP(:), IPX(:), IPY(:), IPZ(:)
      Integer ITYPE, NIP, IC_End, IC_Str, IC, JOPT, ICMP, IRC, I, ISCHK

C NOW DO INTEGRATION WITH AO MATRICES
C FOR THE EXPECTATION VALUE

C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

C ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
c The extra 4 elements correspond to the nuclear contribution
c and the origin of the operator
      NIP=4+NBTRI
      CALL mma_allocate(IP,NIP,Label='IP')
      If (iOpt.eq.1) Then
         CALL mma_allocate(IPX,NIP,Label='IPX')
         CALL mma_allocate(IPY,NIP,Label='IPY')
         CALL mma_allocate(IPZ,NIP,Label='IPZ')
         IPX(:)=0.0D0
         IPY(:)=0.0D0
         IPZ(:)=0.0D0
      End If

      If (iOpt.eq.1) Then
         IC_End=3
         IC_Str=1
      Else
         IC_End=IC_
         IC_Str=IC_
      End If
      DO IC=IC_Str,IC_End ! loop over reading X,Y, and Z AO Integrals

C Get info from the stored integrals
c JOPT controls what is read.
c JOPT=1 Read the size information
c JOPT=0 Read the property
c JOPT=6 Read the property, skipping the nuclear contribution and the origin
c (see OneDat module)
      JOPT=ibset(0,sOpSiz)
      ICMP=IC
      CALL iRDONE(IRC,JOPT,CHARPROP,ICMP,IDUM,ISCHK)

c Actually read the integral
      JOPT=0
      CALL RDONE(IRC,JOPT,CHARPROP,ICMP,IP,ISCHK)

      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)')'*** ERROR IN SUBROUTINE SONATORB ***'
        WRITE(6,'(6X,A)')'  FAILED IN READING FROM  ONEINT'
        WRITE(6,'(6X,A,A)')'  LABEL     = ',CHARPROP
        WRITE(6,'(6X,A,I2)')'  COMPONENT = ',IC
        WRITE(6,*)
        CALL ABEND()
      END IF

      If (iOpt.eq.1) Then
c        note reordering
             CALL DAXPY_(NIP,ROTMAT(IC,1),IP,1,IPX,1)
             CALL DAXPY_(NIP,ROTMAT(IC,2),IP,1,IPY,1)
             CALL DAXPY_(NIP,ROTMAT(IC,3),IP,1,IPZ,1)
      End If

      END DO ! end loop over reading X,Y, and Z AO Integrals

      PROPVALXR=0.0d0
      PROPVALYR=0.0d0
      PROPVALZR=0.0d0
      PROPVALXI=0.0d0
      PROPVALYI=0.0d0
      PROPVALZI=0.0d0

C The integral is NBTRI matrix
C The property is NBTRI matrix
c We only work with half the matrix. Therefore, this would
c have a factor of 2. However, the factor of 1/2 was missing
c in SONATORB.F from the symmetric/antisymmetric equations
      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        If (iOpt.eq.1) Then
        DO I=1,NBTRI
          PROPVALXR=PROPVALXR+IPX(I)*DENS(1,I)
          PROPVALYR=PROPVALYR+IPY(I)*DENS(2,I)
          PROPVALZR=PROPVALZR+IPZ(I)*DENS(3,I)

          PROPVALXI=PROPVALXI+IPX(I)*DENS(4,I)
          PROPVALYI=PROPVALYI+IPY(I)*DENS(5,I)
          PROPVALZI=PROPVALZI+IPZ(I)*DENS(6,I)
        END DO
        Else
        DO I=1,NBTRI
          PROPVALXR=PROPVALXR+IP(I)*DENS(1,I)
          PROPVALYR=PROPVALYR+IP(I)*DENS(2,I)
          PROPVALZR=PROPVALZR+IP(I)*DENS(3,I)

          PROPVALXI=PROPVALXI+IP(I)*DENS(4,I)
          PROPVALYI=PROPVALYI+IP(I)*DENS(5,I)
          PROPVALZI=PROPVALZI+IP(I)*DENS(6,I)
        END DO
        End If
      ELSE
        If (iOpt.eq.1) Then
        DO I=1,NBTRI
          PROPVALXI=PROPVALXI+IPX(I)*DENS(1,I)
          PROPVALYI=PROPVALYI+IPY(I)*DENS(2,I)
          PROPVALZI=PROPVALZI+IPZ(I)*DENS(3,I)

          PROPVALXR=PROPVALXR-IPX(I)*DENS(4,I)
          PROPVALYR=PROPVALYR-IPY(I)*DENS(5,I)
          PROPVALZR=PROPVALZR-IPZ(I)*DENS(6,I)
        END DO
        Else
        DO I=1,NBTRI
          PROPVALXI=PROPVALXI+IP(I)*DENS(1,I)
          PROPVALYI=PROPVALYI+IP(I)*DENS(2,I)
          PROPVALZI=PROPVALZI+IP(I)*DENS(3,I)

          PROPVALXR=PROPVALXR-IP(I)*DENS(4,I)
          PROPVALYR=PROPVALYR-IP(I)*DENS(5,I)
          PROPVALZR=PROPVALZR-IP(I)*DENS(6,I)
        END DO
        End If
      END IF

      WRITE(6,*)
      WRITE(6,*) "************************************"
      WRITE(6,*) "SONATORB EXPECTATION VALUES"
      WRITE(6,*) " PROPERTY: ", CHARPROP
      WRITE(6,*) " TYPE: ", CHARTYPE
      WRITE(6,*) " STATE (K,L): ",ASS,BSS
      WRITE(6,*) "************************************"
      WRITE(6,*) "Property: Re(X): ",PROPVALXR
      WRITE(6,*) "Property: Re(Y): ",PROPVALYR
      WRITE(6,*) "Property: Re(Z): ",PROPVALZR
      WRITE(6,*) "Property: Im(X): ",PROPVALXI
      WRITE(6,*) "Property: Im(Y): ",PROPVALYI
      WRITE(6,*) "Property: Im(Z): ",PROPVALZI
      WRITE(6,*) "************************************"

c Free up un-needed space
      call mma_deallocate(IP)
      If (iOpt.eq.1) Then
         call mma_deallocate(IPX)
         call mma_deallocate(IPY)
         call mma_deallocate(IPZ)
      End If

      END SUBROUTINE SONATORBM_INT

