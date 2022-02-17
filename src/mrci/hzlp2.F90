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

use mrci_global, only: ESMALL, IDFREE, IDISKC, IDISKD, IDISKS, LUEIG, MBUF, MXVEC, MXZ, NCONF, NNEW, NRROOT, NSECT, NVEC, NVTOT, &
                       VSMALL, VZERO
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: CBUF(MBUF,MXVEC), SBUF(MBUF,MXVEC), DBUF(MBUF), CSECT(NSECT,MXVEC), RSECT(NSECT,MXVEC), &
                              XI1(NSECT,NRROOT), XI2(NSECT,NRROOT), CNEW(NSECT,NRROOT)
integer(kind=iwp), intent(out) :: ICI(MBUF)
integer(kind=iwp) :: I, IBUF, IDD, IEND, ISECT, ISTA, IVZ1, IVZ2, IVZ3, IVZ4, IVZSTA, JEND, JSTA, K, KK, NN, NNVEC
integer(kind=iwp), allocatable :: IDCR(:), IDCW(:), IDS(:)

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
!write(u6,*)
!write(u6,*) ' IN HZLP2. NNEW=',NNEW
!if (NVEC < MXVEC) write(u6,*) ' DUMMY WRITES OF NEW FUNCTIONS.'
! WE MAY NEED DUMMY WRITES TO PROVIDE DISK ADDRESSES:
NNVEC = min(NVEC+NNEW,MXVEC)
do K=NVEC+1,NNVEC
  IDISKC(K) = IDFREE
  !write(u6,'(A,I2,A,I8)') ' IDISKC(',K,')=',IDFREE
  do ISTA=1,NCONF,MBUF
    IEND = min(NCONF,ISTA+MBUF-1)
    IBUF = 1+IEND-ISTA
    call iDAFILE(LUEIG,0,ICI,IBUF,IDFREE)
  end do
end do
! WE NEED COPIES OF THE DISK ADDRESSES. TWO COPIES FOR PSI BUFFERS.
call mma_allocate(IDCR,NNVEC,label='IDCR')
call mma_allocate(IDCW,NNVEC,label='IDCW')
call mma_allocate(IDS,NNVEC,label='IDS')
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
    call DGEMM_('N','N',ISECT,NRROOT,NVEC,One,CBUF(JSTA,1),MBUF,VSMALL,MXVEC,Zero,CSECT,NSECT)
    ! THEN, SIGMA SECTION INTO RSECT.
    call DGEMM_('N','N',ISECT,NRROOT,NVEC,One,SBUF(JSTA,1),MBUF,VSMALL,MXVEC,Zero,RSECT,NSECT)
    ! AND THEN FORM RSECT=SECTION OF RESIDUAL ARRAY, AND XI1 AND XI2:
    do I=1,ISECT
      do K=1,NRROOT
        RSECT(I,K) = RSECT(I,K)-ESMALL(K)*CSECT(I,K)
        XI1(I,K) = CSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
        XI2(I,K) = RSECT(I,K)/(DBUF(I+JSTA-1)-ESMALL(K))
      end do
    end do
    ! FORM NEW PSI ARRAYS IN CNEW SECTION:
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,One,CSECT,NSECT,VZERO(IVZ1,IVZSTA),MXZ,Zero,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,One,RSECT,NSECT,VZERO(IVZ2,IVZSTA),MXZ,One,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,One,XI1,NSECT,VZERO(IVZ3,IVZSTA),MXZ,One,CNEW,NSECT)
    call DGEMM_('N','N',ISECT,NNEW,NRROOT,One,XI2,NSECT,VZERO(IVZ4,IVZSTA),MXZ,One,CNEW,NSECT)
    !if (ISTA+JSTA == 2) then
    !  write(u6,*) ' CONSTRUCTION OF NEW VECTOR IN HZLP2.'
    !  write(u6,*) ' CSECT:'
    !  write(u6,'(1X,5F15.6)') ((CSECT(I,J),I=1,5),J=1,NNEW)
    !  write(u6,*) ' RSECT:'
    !  write(u6,'(1X,5F15.6)') ((RSECT(I,J),I=1,5),J=1,NNEW)
    !  write(u6,*) '   XI1:'
    !  write(u6,'(1X,5F15.6)') ((XI1(I,J),I=1,5),J=1,NNEW)
    !  write(u6,*) '   XI2:'
    !  write(u6,'(1X,5F15.6)') ((XI2(I,J),I=1,5),J=1,NNEW)
    !  write(u6,*) ' VZERO:'
    !  IIII = IVZSTA-1
    !  write(u6,'(1X,4F15.6)') ((VZERO(I,IIII+J),I=1,4*NNEW),J=1,NNEW)
    !  write(u6,*) '  CNEW:'
    !  write(u6,'(1X,5F15.6)') ((CNEW(I,J),I=1,5),J=1,NNEW)
    !end if
    ! INSERT THE NEW PSI SECTIONS IN BUFFER. THIS MAY IMPLY OVERWRITING
    ! OLD ENTRIES, BUT CAN ALSO LEAD TO AN INCREASED NUMBER OF VECTORS:
    do KK=1,NNEW
      NN = NVTOT+KK
      K = 1+mod(NN-1,MXVEC)
      !if (ISTA+JSTA == 2) then
      !  write(u6,'(A,I2,A,I6)') ' CNEW NR.',KK,' COPIED TO BUFFER ',K
      !  write(u6,*) ' IT CONTAINS:'
      !  write(u6,'(1X,5F15.6)') (CNEW(I,KK),I=1,15)
      !end if
      CBUF(JSTA:JSTA+ISECT-1,K) = CNEW(1:ISECT,KK)
    end do
    ! CONTINUE, NEXT SECTION.
  end do
  do KK=1,NNEW
    NN = NVTOT+KK
    K = 1+mod(NN-1,MXVEC)
    !if (ISTA == 1) then
    !  write(u6,'(A,I2,A,I6)') ' BUFFER NR.',K,' WRITTEN AT ',IDCW(K)
    !  write(u6,*) ' IT CONTAINS:'
    !  write(u6,'(1X,5F15.6)') (CBUF(I,K),I=1,15)
    !end if
    call PKVEC(IBUF,CBUF(1,K),ICI)
    call iDAFILE(LUEIG,1,ICI,IBUF,IDCW(K))
  end do
  ! CONTINUE, NEXT BUFFER.
end do
NVTOT = NVTOT+NNEW
call mma_deallocate(IDCR)
call mma_deallocate(IDCW)
call mma_deallocate(IDS)

return

end subroutine HZLP2
