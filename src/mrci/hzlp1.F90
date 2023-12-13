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

subroutine HZLP1(CBUF,SBUF,DBUF,ARR,CSECT,RSECT,XI1,XI2,ICI)

use mrci_global, only: ESMALL, IDISKC, IDISKD, IDISKS, LUEIG, MBUF, MXVEC, NCONF, NRROOT, NSECT, NVEC, VSMALL
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: CBUF(MBUF,MXVEC), SBUF(MBUF,MXVEC), DBUF(MBUF), ARR(NRROOT,NRROOT,11), CSECT(NSECT,MXVEC), &
                              RSECT(NSECT,MXVEC), XI1(NSECT,NRROOT), XI2(NSECT,NRROOT)
integer(kind=iwp), intent(out) :: ICI(MBUF)
integer(kind=iwp) :: I, IBUF, IDD, IEND, ISECT, ISTA, JEND, JSTA, K
integer(kind=iwp), parameter :: IX1F = 1, IX2F = 2, IRR = 3, IX1R = 4, IX2R = 5, IX1X1 = 6, IX2X1 = 7, IX2X2 = 8, IFDF = 9, &
                                IFDR = 10, IRDR = 11
integer(kind=iwp), allocatable :: IDC(:), IDS(:)

! THIS SUBROUTINE LOOPS OVER SECTIONS OF PSI AND SIGMA ARRAYS
! ON DISK, AND ACCUMULATES OVERLAP MATRICES AND A COUPLE OF
! HAMILTONIAN MATRICES IN THE BASIS SET PSI, RHO, XI1 AND XI2.
! THE 11 MATRICES X1F,..,RDR ARE STORED CONSECUTIVELY IN THE
! SINGLE ARRAY ARR.
ARR(:,:,:) = Zero
call mma_allocate(IDC,NVEC,label='IDC')
call mma_allocate(IDS,NVEC,label='IDS')
do K=1,NVEC
  IDC(K) = IDISKC(K)
  IDS(K) = IDISKS(K)
end do
IDD = IDISKD
! LOOP OVER BUFFERS FOR READING PSI, SIGMA AND DBUF:
do ISTA=1,NCONF,MBUF
  IEND = min(NCONF,ISTA+MBUF-1)
  IBUF = 1+IEND-ISTA
  call dDAFILE(LUEIG,2,DBUF,IBUF,IDD)
  do K=1,NVEC
    call iDAFILE(LUEIG,2,ICI,IBUF,IDC(K))
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
    ! ACCUMULATE OVERLAP MATRICES:
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI1,NSECT,CSECT,NSECT,One,ARR(1,1,IX1F),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI2,NSECT,CSECT,NSECT,One,ARR(1,1,IX2F),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,RSECT,NSECT,RSECT,NSECT,One,ARR(1,1,IRR),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI1,NSECT,RSECT,NSECT,One,ARR(1,1,IX1R),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI2,NSECT,RSECT,NSECT,One,ARR(1,1,IX2R),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI1,NSECT,XI1,NSECT,One,ARR(1,1,IX1X1),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI2,NSECT,XI1,NSECT,One,ARR(1,1,IX2X1),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI2,NSECT,XI2,NSECT,One,ARR(1,1,IX2X2),NRROOT)
    ! PUT D*CSECT INTO XI1, AND D*RSECT INTO XI2:
    do I=1,ISECT
      do K=1,NRROOT
        XI1(I,K) = DBUF(I+JSTA-1)*CSECT(I,K)
        XI2(I,K) = DBUF(I+JSTA-1)*RSECT(I,K)
      end do
    end do
    ! ACCUMULATE ARRAYS FDF, FDR, AND RDR:
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI1,NSECT,CSECT,NSECT,One,ARR(1,1,IFDF),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI1,NSECT,RSECT,NSECT,One,ARR(1,1,IFDR),NRROOT)
    call DGEMM_('T','N',NRROOT,NRROOT,ISECT,One,XI2,NSECT,RSECT,NSECT,One,ARR(1,1,IRDR),NRROOT)
    ! CONTINUE, NEXT SECTION.
  end do
  ! CONTINUE, NEXT BUFFER.
end do
call mma_deallocate(IDC)
call mma_deallocate(IDS)

return

end subroutine HZLP1
