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

subroutine HZLP2(CBUF,SBUF,DBUF,CSECT,RSECT,XI1,XI2,CNEW,ICI)

implicit real*8(A-H,O-Z)
parameter(ONE=1.0d00)
#include "SysDef.fh"
#include "mrci.fh"
dimension CBUF(MBUF,MXVEC), SBUF(MBUF,MXVEC), DBUF(MBUF), ICI(MBUF)
dimension CSECT(NSECT,MXVEC), RSECT(NSECT,MXVEC)
dimension CNEW(NSECT,MXVEC)
dimension XI1(NSECT,MXVEC), XI2(NSECT,MXVEC)
dimension IDCR(MXVEC), IDCW(MXVEC), IDS(MXVEC)

! THIS SUBROUTINE LOOPS OVER SECTIONS OF PSI AND SIGMA ARRAYS
! ON DISK, AND FORMS A NEW SET OF PSI ARRAYS AS A LINEAR COMBINATION
! OF THE BASIS SET PSI, RHO, XI1 AND XI2. TO FORM THE NEW PSI
! ARRAY, THE FIRST (NRROOT-NNEW) COLUMNS ARE SKIPPED, SINCE THEY
! GIVE NO ESSENTIAL IMPROVEMENT.
IVZSTA = 1+NRROOT-NNEW
IVZ1 = 1
IVZ2 = 1+NRROOT
IVZ3 = 1+2*NRROOT
IVZ4 = 1+3*NRROOT
!write(6,*)
!write(6,*) ' IN HZLP2. NNEW=',NNEW
!if (NVEC < MXVEC) write(6,*) ' DUMMY WRITES OF NEW FUNCTIONS.'
! WE MAY NEED DUMMY WRITES TO PROVIDE DISK ADDRESSES:
NNVEC = min(NVEC+NNEW,MXVEC)
do K=NVEC+1,NNVEC
  IDISKC(K) = IDFREE
  !write(6,'(A,I2,A,I8)') ' IDISKC(',K,')=',IDFREE
  do ISTA=1,NCONF,MBUF
    IEND = min(NCONF,ISTA+MBUF-1)
    IBUF = 1+IEND-ISTA
    call iDAFILE(LUEIG,0,ICI,IBUF,IDFREE)
  end do
end do
! WE NEED COPIES OF THE DISK ADDRESSES. TWO COPIES FOR PSI BUFFERS.
do K=1,NNVEC
  IDCR(K) = IDISKC(K)
  IDCW(K) = IDISKC(K)
  IDS(K) = IDISKS(K)
end do
IDD = IDISKD
! LOOP OVER BUFFERS FOR READING PSI, SIGMA AND DBUF:
do ISTA=1,NCONF,MBUF
  IEND = min(NCONF,ISTA+MBUF-1)
  IBUF = 1+IEND-ISTA
  call dDAFILE(LUEIG,2,DBUF,IBUF,IDD)
  do K=1,NVEC
    call iDAFILE(LUEIG,2,ICI,IBUF,IDCR(K))
    call UPKVEC(IBUF,ICI,CBUF(1,K))
    call dDAFILE(LUEIG,2,SBUF(1,K),IBUF,IDS(K))
  end do
  ! LOOP OVER VECTOR SECTIONS, LENGTH AT MOST NSECT:
  do JSTA=1,IBUF,NSECT
    JEND = min(IBUF,JSTA+NSECT-1)
    ISECT = 1+JEND-JSTA
    ! TRANSFORM TO EIGENFUNCTIONS OF HSMALL: FIRST, CI SECTION.
    call DGEMM_('N','N',ISECT,NRROOT,NVEC,1.0d0,CBUF(JSTA,1),MBUF,VSMALL,MXVEC,0.0d0,CSECT,NSECT)
    ! THEN, SIGMA SECTION INTO RSECT.
    call DGEMM_('N','N',ISECT,NRROOT,NVEC,1.0d0,SBUF(JSTA,1),MBUF,VSMALL,MXVEC,0.0d0,RSECT,NSECT)
    ! AND THEN FORM RSECT=SECTION OF RESIDUAL ARRAY, AND XI1 AND XI2:
    do I=1,ISECT
      do K=1,NRROOT
        RSECT(I,K) = RSECT(I,K)-ESMALL(K)*CSECT(I,K)
        XI1(I,K) = CSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
        XI2(I,K) = RSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
      end do
    end do
    ! FORM NEW PSI ARRAYS IN CNEW SECTION:
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,1.0d0,CSECT,NSECT,VZERO(IVZ1,IVZSTA),MXZ,0.0d0,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,RSECT,NSECT,VZERO(IVZ2,IVZSTA),MXZ,ONE,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,XI1,NSECT,VZERO(IVZ3,IVZSTA),MXZ,ONE,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,ONE,XI2,NSECT,VZERO(IVZ4,IVZSTA),MXZ,ONE,CNEW,NSECT)
    !if (ISTA+JSTA == 2) then
    !  write(6,*) ' CONSTRUCTION OF NEW VECTOR IN HZLP2.'
    !  write(6,*) ' CSECT:'
    !  write(6,'(1X,5F15.6)') ((CSECT(I,J),I=1,5),J=1,NNEW)
    !  write(6,*) ' RSECT:'
    !  write(6,'(1X,5F15.6)') ((RSECT(I,J),I=1,5),J=1,NNEW)
    !  write(6,*) '   XI1:'
    !  write(6,'(1X,5F15.6)') ((XI1(I,J),I=1,5),J=1,NNEW)
    !  write(6,*) '   XI2:'
    !  write(6,'(1X,5F15.6)') ((XI2(I,J),I=1,5),J=1,NNEW)
    !  write(6,*) ' VZERO:'
    !  IIII = IVZSTA-1
    !  write(6,'(1X,4F15.6)') ((VZERO(I,IIII+J),I=1,4*NNEW),J=1,NNEW)
    !  write(6,*) '  CNEW:'
    !  write(6,'(1X,5F15.6)') ((CNEW(I,J),I=1,5),J=1,NNEW)
    !end if
    ! INSERT THE NEW PSI SECTIONS IN BUFFER. THIS MAY IMPLY OVERWRITING
    ! OLD ENTRIES, BUT CAN ALSO LEAD TO AN INCREASED NUMBER OF VECTORS:
    do KK=1,NNEW
      NN = NVTOT+KK
      K = 1+mod(NN-1,MXVEC)
      !if (ISTA+JSTA == 2) then
      !  write(6,'(A,I2,A,I6)') ' CNEW NR.',KK,' COPIED TO BUFFER ',K
      !  write(6,*) ' IT CONTAINS:'
      !  write(6,'(1X,5F15.6)') (CNEW(I,KK),I=1,15)
      !end if
      call DCOPY_(ISECT,CNEW(1,KK),1,CBUF(JSTA,K),1)
    end do
  ! CONTINUE, NEXT SECTION.
  end do
  do KK=1,NNEW
    NN = NVTOT+KK
    K = 1+mod(NN-1,MXVEC)
    !if (ISTA == 1) then
    !  write(6,'(A,I2,A,I6)') ' BUFFER NR.',K,' WRITTEN AT ',IDCW(K)
    !  write(6,*) ' IT CONTAINS:'
    !  write(6,'(1X,5F15.6)') (CBUF(I,K),I=1,15)
    !end if
    call PKVEC(IBUF,CBUF(1,K),ICI)
    call iDAFILE(LUEIG,1,ICI,IBUF,IDCW(K))
  end do
  ! CONTINUE, NEXT BUFFER.
end do
NVTOT = NVTOT+NNEW

return

end subroutine HZLP2
