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
      SUBROUTINE DENS(IVEC,DMAT)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"
#include "para_info.fh"
      DIMENSION DMAT(*)

      CALL QENTER('DENS')
C Compute total density matrix as symmetry-blocked array of
C triangular matrices in DMAT. Size of a triangular submatrix is
C  (NORB(ISYM)*(NORB(ISYM)+1))/2.
      NDMAT=0
      NDPT=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        NDPT=NDPT+NO**2
        NDMAT=NDMAT+(NO*(NO+1))/2
      END DO
      CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
C First, put in the reference density matrix.
      IDMOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO II=1,NI
          IDM=IDMOFF+(II*(II+1))/2
          DMAT(IDM)=2.0D0
        END DO
        DO IT=1,NA
          ITABS=NAES(ISYM)+IT
          ITTOT=NI+IT
          DO IU=1,IT
            IUABS=NAES(ISYM)+IU
            IUTOT=NI+IU
            IDRF=(ITABS*(ITABS-1))/2+IUABS
            IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
            DMAT(IDM)=WORK(LDREF-1+IDRF)
          END DO
        END DO
         IDMOFF=IDMOFF+(NO*(NO+1))/2
      END DO
*     WRITE(*,*)' DENS. Initial DMAT:'
*     WRITE(*,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
      CALL GETMEM('DPT','ALLO','REAL',LDPT,NDPT)
      CALL GETMEM('DSUM','ALLO','REAL',LDSUM,NDPT)
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDSUM),1)

C The 1st order contribution to the density matrix
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
      CALL TRDNS1(IVEC,WORK(LDPT))
      CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*     WRITE(*,*)' DPT after TRDNS1.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
      CALL TRDNS2D(IVEC,IVEC,WORK(LDPT),NDPT)
      IF(IFDENS) THEN
C The exact density matrix evaluation:
        CALL TRDTMP(WORK(LDPT))
      ELSE
C The approximate density matrix evaluation:
        CALL TRDNS2A(IVEC,IVEC,WORK(LDPT))
      END IF
      CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*     WRITE(*,*)' DPT after TRDNS2D.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
      CALL TRDNS2O(IVEC,IVEC,WORK(LDPT))
      CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*     WRITE(*,*)' DPT after TRDNS2O.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      CALL GETMEM('DPT','FREE','REAL',LDPT,NDPT)
      IDMOFF=0
      IDSOFF=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        DO IP=1,NO
          DO IQ=1,IP
            IDM=IDMOFF+(IP*(IP-1))/2+IQ
            IDSUM=IDSOFF+IP+NO*(IQ-1)
            DMAT(IDM)=DMAT(IDM)+WORK(LDSUM-1+IDSUM)
          END DO
        END DO
        IDMOFF=IDMOFF+(NO*(NO+1))/2
        IDSOFF=IDSOFF+NO**2
      END DO
      CALL GETMEM('DSUM','FREE','REAL',LDSUM,NDPT)
C Scale with 1/DENORM to normalize
      X=1.0D0/DENORM
      CALL DSCAL_(NDMAT,X,DMAT,1)

CSVC: For true parallel calculations, replicate the DMAT array
C so that the slaves have the same density matrix as the master.
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        IF (.NOT.KING()) THEN
          CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
        END IF
        CALL GADSUM(DMAT,NDMAT)
      END IF
#endif

      CALL QEXIT('DENS')
      RETURN
      END
