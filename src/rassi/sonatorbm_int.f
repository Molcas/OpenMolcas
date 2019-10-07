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
      SUBROUTINE SONATORBM_INT(DENS, CHARPROP, IC_,CHARTYPE,ASS,BSS,NSS,
     &                         iOpt,ROTMAT,
     &                         PROPVALXR,PROPVALYR,PROPVALZR,
     &                         PROPVALXI,PROPVALYI,PROPVALZI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORBM_INT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      INTEGER ASS,BSS
      CHARACTER*8 CHARPROP, CHARTYPE
      Dimension ROTMAT(3,3)
      DIMENSION IDUM(1)

C NOW DO INTEGRATION WITH AO MATRICES
C FOR THE EXPECTATION VALUE

C The following creates an array that is used to
C map a specific spin state to the corresponding
C spin-free state and to its spin
C (see prprop.f and others)

      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)

      ISS=0
      DO ISF=1,NSTATE
        JOB=iWork(lJBNUM+ISF-1)
        MPLET=MLTPLT(JOB)

        DO MSPROJ=-MPLET+1,MPLET-1,2
          ISS=ISS+1
          IWORK(LMAPST-1+ISS)=ISF
          IWORK(LMAPSP-1+ISS)=MPLET
          IWORK(LMAPMS-1+ISS)=MSPROJ
        END DO
      END DO


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
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
      If (iOpt.eq.1) Then
         CALL GETMEM('IPX    ','ALLO','REAL',LIPX,NIP)
         CALL GETMEM('IPY    ','ALLO','REAL',LIPY,NIP)
         CALL GETMEM('IPZ    ','ALLO','REAL',LIPZ,NIP)
         CALL DCOPY_(NIP,[0.0D00],0,WORK(LIPX),1)
         CALL DCOPY_(NIP,[0.0D00],0,WORK(LIPY),1)
         CALL DCOPY_(NIP,[0.0D00],0,WORK(LIPZ),1)
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
c (see misc_util/OneFlags.fh)
      JOPT=1
      CALL iRDONE(IRC,JOPT,CHARPROP,IC,IDUM,ISCHK)
      IF (IRC.EQ.0) NSIZ=IDUM(1)


c Actually read the integral
      JOPT=0
      CALL RDONE(IRC,JOPT,CHARPROP,IC,WORK(LIP),ISCHK)

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
             CALL DAXPY_(NIP,ROTMAT(IC,1),WORK(LIP),1,WORK(LIPX),1)
             CALL DAXPY_(NIP,ROTMAT(IC,2),WORK(LIP),1,WORK(LIPY),1)
             CALL DAXPY_(NIP,ROTMAT(IC,3),WORK(LIP),1,WORK(LIPZ),1)
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
          PROPVALXR=PROPVALXR+WORK(LIPX-1+I)*DENS(1,I)
          PROPVALYR=PROPVALYR+WORK(LIPY-1+I)*DENS(2,I)
          PROPVALZR=PROPVALZR+WORK(LIPZ-1+I)*DENS(3,I)

          PROPVALXI=PROPVALXI+WORK(LIPX-1+I)*DENS(4,I)
          PROPVALYI=PROPVALYI+WORK(LIPY-1+I)*DENS(5,I)
          PROPVALZI=PROPVALZI+WORK(LIPZ-1+I)*DENS(6,I)
        END DO
        Else
        DO I=1,NBTRI
          PROPVALXR=PROPVALXR+WORK(LIP-1+I)*DENS(1,I)
          PROPVALYR=PROPVALYR+WORK(LIP-1+I)*DENS(2,I)
          PROPVALZR=PROPVALZR+WORK(LIP-1+I)*DENS(3,I)

          PROPVALXI=PROPVALXI+WORK(LIP-1+I)*DENS(4,I)
          PROPVALYI=PROPVALYI+WORK(LIP-1+I)*DENS(5,I)
          PROPVALZI=PROPVALZI+WORK(LIP-1+I)*DENS(6,I)
        END DO
        End If
      ELSE
        If (iOpt.eq.1) Then
        DO I=1,NBTRI
          PROPVALXI=PROPVALXI+WORK(LIPX-1+I)*DENS(1,I)
          PROPVALYI=PROPVALYI+WORK(LIPY-1+I)*DENS(2,I)
          PROPVALZI=PROPVALZI+WORK(LIPZ-1+I)*DENS(3,I)

          PROPVALXR=PROPVALXR-WORK(LIPX-1+I)*DENS(4,I)
          PROPVALYR=PROPVALYR-WORK(LIPY-1+I)*DENS(5,I)
          PROPVALZR=PROPVALZR-WORK(LIPZ-1+I)*DENS(6,I)
        END DO
        Else
        DO I=1,NBTRI
          PROPVALXI=PROPVALXI+WORK(LIP-1+I)*DENS(1,I)
          PROPVALYI=PROPVALYI+WORK(LIP-1+I)*DENS(2,I)
          PROPVALZI=PROPVALZI+WORK(LIP-1+I)*DENS(3,I)

          PROPVALXR=PROPVALXR-WORK(LIP-1+I)*DENS(4,I)
          PROPVALYR=PROPVALYR-WORK(LIP-1+I)*DENS(5,I)
          PROPVALZR=PROPVALZR-WORK(LIP-1+I)*DENS(6,I)
        END DO
        End If
      END IF

      IF(IPGLOB.GE.VERBOSE) THEN
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
      END IF

c Free up un-needed space
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      If (iOpt.eq.1) Then
         CALL GETMEM('IPX    ','FREE','REAL',LIPX,NIP)
         CALL GETMEM('IPY    ','FREE','REAL',LIPY,NIP)
         CALL GETMEM('IPZ    ','FREE','REAL',LIPZ,NIP)
      End If

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      RETURN
      END

