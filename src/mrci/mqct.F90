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

subroutine MQCT(AREF,EREF,CI,SGM,ICI)

use mrci_global, only: CSPCK, ENGY, ESHIFT, ESMALL, ETHRE, FOCK, GFAC, HZERO, ICPF, IDFREE, IDISKC, IDISKS, INDX, INTSY, IPRINT, &
                       IREFX, IROOT, ISAB, ISMAX, ITER, JREFX, KBUFF1, LUEIG, LUREST, MAXIT, MBUF, MXVEC, MXZ, NBMN, NCONF, NNEW, &
                       NREF, NRROOT, NSECT, NSTOT, NVEC, NVMAX, NVSQ, NVTOT, SQNLIM, SZERO, VSMALL, VZERO
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: AREF(NREF,NREF), EREF(NREF)
real(kind=wp), intent(out) :: CI(NCONF), SGM(NCONF)
integer(kind=iwp), intent(out) :: ICI(MBUF)
#include "warnings.h"
integer(kind=iwp) :: I, IBUF, ICSF, IDISK, IDREST, IEND, II, IMAX, IMIN, IPOS, IR, IRR, ISTA, IVEC, J, K, KK, KL, L, LL, NCONV, &
                     NN, NRON, NZ
real(kind=wp) :: C, C2NREF, C2REF, CPTIT, CPTNOW, CPTOLD, CPTOT, CPTSTA, DUM1, DUM2, DUM3, EACPF, ECI, EDAV, EDISP, ELOW, EMIN, &
                 ENREF, H, P, PMAX, QACPF, QDAV, RSUM, S, SQNRM, THR, TMP
integer(kind=iwp), allocatable :: IBMN(:), IDC(:), IDS(:)
real(kind=wp), allocatable :: ABIJ(:), AC1(:), AC2(:), AIBJ(:), AJBI(:), ARR(:,:,:), ASCR1(:), ASCR2(:), BFIN3(:), BFIN4(:), &
                              BIAC2(:), BICA2(:), BMN(:), BSCR1(:), BSCR2(:), CBUF(:,:), CNEW(:,:), CSECT(:,:), DBK(:), DBUF(:), &
                              ELAST(:), EZERO(:), FSCR1(:), FSCR2(:), FSEC(:), HCOPY(:,:), HSMALL(:,:), PCOPY(:,:), PSEL(:), &
                              PSMALL(:,:), RNRM(:), RSECT(:,:), SBUF(:,:), SCOPY(:,:), SCR(:), SSMALL(:,:), XI1(:,:), XI2(:,:)
real(kind=wp), external :: DDOT_

call mma_allocate(CBUF,MBUF,MXVEC,label='CBUF')
call mma_allocate(SBUF,MBUF,MXVEC,label='SBUF')
call mma_allocate(DBUF,MBUF,label='DBUF')
call mma_allocate(CSECT,NSECT,MXVEC,label='CSECT')
call mma_allocate(RSECT,NSECT,MXVEC,label='RSECT')
call mma_allocate(XI1,NSECT,NRROOT,label='XI1')
call mma_allocate(XI2,NSECT,NRROOT,label='XI2')
call mma_allocate(CNEW,NSECT,NRROOT,label='CNEW')
call mma_allocate(SCR,MXVEC*MBUF,label='SCR')
call mma_allocate(IDC,MXVEC,label='IDC')
call mma_allocate(IDS,MXVEC,label='IDS')
call mma_allocate(ELAST,NRROOT,label='ELAST')
call mma_allocate(HCOPY,MXVEC,MXVEC,label='HCOPY')
call mma_allocate(PCOPY,MXVEC,MXVEC,label='PCOPY')
call mma_allocate(SCOPY,MXVEC,MXVEC,label='SCOPY')
call mma_allocate(PSEL,MXVEC,label='PSEL')
call mma_allocate(RNRM,NRROOT,label='RNRM')
call mma_allocate(HSMALL,MXVEC,MXVEC,label='HSMALL')
call mma_allocate(SSMALL,MXVEC,MXVEC,label='SSMALL')
call mma_allocate(PSMALL,MXVEC,MXVEC,label='PSMALL')
call mma_allocate(EZERO,MXZ,label='EZERO')
call mma_allocate(HZERO,MXZ,MXZ,label='HZERO')
call mma_allocate(SZERO,MXZ,MXZ,label='SZERO')
call mma_allocate(VZERO,MXZ,MXZ,label='VZERO')

write(u6,*)
write(u6,*) repeat('-',60)
if (ICPF == 0) then
  write(u6,*) '   MR SDCI CALCULATION.'
else
  write(u6,*) '   MR ACPF CALCULATION.'
end if
write(u6,*) repeat('-',60)
write(u6,*)
write(u6,*) '         CONVERGENCE STATISTICS:'
write(u6,'(1X,A)') 'ITER NVEC     ENERGIES    LOWERING RESIDUAL SEL.WGT CPU(S) CPU TOT'
ITER = 0
call SETTIM()
call TIMING(CPTNOW,DUM1,DUM2,DUM3)
CPTOLD = CPTNOW
CPTSTA = CPTNOW
HSMALL(1,1) = Zero
! LOOP HEAD FOR CI ITERATIONS:
do
  ITER = ITER+1
  ! --------------------------------------------------------------------
  ! CALCULATE SIGMA ARRAYS FOR SHIFTED HAMILTONIAN IN MCSF BASIS:
  do I=1,NNEW
    IVEC = 1+mod(NVTOT-NNEW+I-1,MXVEC)
    IDISK = IDISKC(IVEC)
    do ISTA=1,NCONF,MBUF
      NN = min(MBUF,NCONF+1-ISTA)
      call iDAFILE(LUEIG,2,ICI,NN,IDISK)
      call UPKVEC(NN,ICI,CI(ISTA))
    end do
    call mma_allocate(BMN,NBMN,label='BMN')
    call mma_allocate(IBMN,NBMN,label='IBMN')
    call mma_allocate(BIAC2,ISMAX,label='BIAC2')
    call mma_allocate(BICA2,ISMAX,label='BICA2')
    call mma_allocate(BFIN3,KBUFF1,label='BFIN3')
    call mma_allocate(AC1,ISMAX,label='AC1')
    call mma_allocate(AC2,ISMAX,label='AC2')
    call mma_allocate(BFIN4,KBUFF1,label='BFIN4')
    call mma_allocate(ABIJ,NVSQ,label='ABIJ')
    call mma_allocate(AIBJ,NVSQ,label='AIBJ')
    call mma_allocate(AJBI,NVSQ,label='AJBI')
    call mma_allocate(ASCR1,NVMAX**2,label='ASCR1')
    call mma_allocate(BSCR1,NVMAX**2,label='BSCR1')
    call mma_allocate(FSCR1,NVSQ,label='FSCR1')
    call mma_allocate(FSEC,2*NVSQ,label='FSEC')
    call mma_allocate(ASCR2,NVMAX**2,label='ASCR2')
    call mma_allocate(BSCR2,NVMAX**2,label='BSCR2')
    call mma_allocate(FSCR2,NVSQ,label='FSCR2')
    call mma_allocate(DBK,2*NVSQ,label='DBK')
    call SIGMA(SGM,AREF,CI,INTSY,INDX,BMN,IBMN,BIAC2,BICA2,BFIN3,ISAB,AC1,AC2,BFIN4,ABIJ,AIBJ,AJBI,ASCR1,BSCR1,FSCR1,FSEC,FOCK, &
               ASCR2,BSCR2,FSCR2,DBK)
    call mma_deallocate(BMN)
    call mma_deallocate(IBMN)
    call mma_deallocate(BIAC2)
    call mma_deallocate(BICA2)
    call mma_deallocate(BFIN3)
    call mma_deallocate(AC1)
    call mma_deallocate(AC2)
    call mma_deallocate(BFIN4)
    call mma_deallocate(ABIJ)
    call mma_deallocate(AIBJ)
    call mma_deallocate(AJBI)
    call mma_deallocate(ASCR1)
    call mma_deallocate(BSCR1)
    call mma_deallocate(FSCR1)
    call mma_deallocate(FSEC)
    call mma_deallocate(ASCR2)
    call mma_deallocate(BSCR2)
    call mma_deallocate(FSCR2)
    call mma_deallocate(DBK)
    NSTOT = NSTOT+1
    ! WRITE IT OUT:
    IVEC = 1+mod(NSTOT-1,MXVEC)
    IDISK = IDISKS(IVEC)
    if (IDISK == -1) IDISK = IDFREE
    do ISTA=1,NCONF,MBUF
      NN = min(MBUF,NCONF+1-ISTA)
      call dDAFILE(LUEIG,1,SGM(ISTA),NN,IDISK)
    end do
    if (IDISK > IDFREE) then
      IDISKS(IVEC) = IDFREE
      IDFREE = IDISK
    end if
  end do
  ! --------------------------------------------------------------------
  ! NR OF VECTORS PRESENTLY RETAINED:
  NVEC = min(MXVEC,NVTOT)
  ! --------------------------------------------------------------------
  ! COPY HSMALL, SSMALL AND PSMALL IN REORDERED FORM, BY AGE:
  do L=NNEW+1,NVEC
    LL = 1+mod(NVTOT-L,MXVEC)
    do K=NNEW+1,NVEC
      KK = 1+mod(NVTOT-K,MXVEC)
      HCOPY(K,L) = HSMALL(KK,LL)
      SCOPY(K,L) = SSMALL(KK,LL)
      PCOPY(K,L) = PSMALL(KK,LL)
    end do
  end do
  ! CLEAR NEW AREAS TO BE USED:
  HCOPY(1:NVEC,1:NNEW) = Zero
  SCOPY(1:NVEC,1:NNEW) = Zero
  PCOPY(1:NVEC,1:NNEW) = Zero
  ! THEN LOOP OVER BUFFERS. FIRST GET COPIES OF DISK ADDRESSES:
  IDC(1:NVEC) = IDISKC(1:NVEC)
  IDS(1:NVEC) = IDISKS(1:NVEC)
  do ISTA=1,NCONF,MBUF
    IEND = min(NCONF,ISTA+MBUF-1)
    IBUF = 1+IEND-ISTA
    do K=1,NVEC
      KK = 1+mod(NVTOT-K,MXVEC)
      call iDAFILE(LUEIG,2,ICI,IBUF,IDC(KK))
      call UPKVEC(IBUF,ICI,CBUF(1,K))
      if (K <= NNEW) call dDAFILE(LUEIG,2,SBUF(1,K),IBUF,IDS(KK))
    end do
    ! ------------------------------------------------------------------
    ! NOTE: AT THIS POINT, THE COLUMNS NR 1..NVEC OF CBUF WILL
    ! CONTAIN THE BUFFERS OF, FIRST, THE NNEW NEWEST PSI ARRAYS,
    ! THEN, THE NOLD ONES FROM EARLIER ITERATIONS.
    ! THE COLUMNS 1..NNEW OF SBUF WILL CONTAIN THE NEWEST NNEW
    ! SIGMA ARRAYS. LEADING DIMENSION OF CBUF AND SBUF IS MBUF. ACTUAL
    ! BUFFER SIZE IS IBUF, WHICH CAN BE SMALLER. ACCUMULATE:
    call DGEMM_('T','N',NVEC,NNEW,IBUF,One,CBUF,MBUF,CBUF,MBUF,Zero,SCR,NVEC)
    KL = 0
    do L=1,NNEW
      do K=1,NVEC
        KL = KL+1
        SCOPY(K,L) = SCOPY(K,L)+SCR(KL)
      end do
    end do
    call DGEMM_('T','N',NVEC,NNEW,IBUF,One,CBUF,MBUF,SBUF,MBUF,Zero,SCR,NVEC)
    KL = 0
    do L=1,NNEW
      do K=1,NVEC
        KL = KL+1
        HCOPY(K,L) = HCOPY(K,L)+SCR(KL)
      end do
    end do
    ! ALSO, UPDATE PSMALL, WHICH IS USED FOR SELECTION.
    if (ISTA > IREFX(NRROOT)) cycle
    do I=1,NRROOT
      IR = IROOT(I)
      IRR = IREFX(IR)
      if ((IRR < ISTA) .or. (IRR > IEND)) cycle
      IPOS = IRR+1-ISTA
      do L=1,NNEW
        PCOPY(1:NVEC,L) = PCOPY(1:NVEC,L)+CBUF(IPOS,1:NVEC)*CBUF(IPOS,L)
      end do
    end do
  end do
  ! TRANSFER ELEMENTS BACK TO HSMALL, ETC.
  do L=1,NNEW
    LL = 1+mod(NVTOT-L,MXVEC)
    do K=1,NVEC
      KK = 1+mod(NVTOT-K,MXVEC)
      H = HCOPY(K,L)
      S = SCOPY(K,L)
      P = PCOPY(K,L)
      HCOPY(L,K) = H
      SCOPY(L,K) = S
      PCOPY(L,K) = P
      HSMALL(KK,LL) = H
      SSMALL(KK,LL) = S
      PSMALL(KK,LL) = P
      HSMALL(LL,KK) = H
      SSMALL(LL,KK) = S
      PSMALL(LL,KK) = P
    end do
  end do
  if (IPRINT >= 10) then
    write(u6,*)
    write(u6,*) ' HSMALL MATRIX:'
    do I=1,NVEC
      write(u6,'(1X,5F15.6)') (HSMALL(I,J),J=1,NVEC)
    end do
    write(u6,*)
    write(u6,*) ' SSMALL MATRIX:'
    do I=1,NVEC
      write(u6,'(1X,5F15.6)') (SSMALL(I,J),J=1,NVEC)
    end do
    write(u6,*)
    write(u6,*) ' PSMALL MATRIX:'
    do I=1,NVEC
      write(u6,'(1X,5F15.6)') (PSMALL(I,J),J=1,NVEC)
    end do
    !write(u6,*)
    !write(u6,*)
    !write(u6,*) ' HCOPY MATRIX:'
    !do I=1,NVEC
    !  write(u6,'(1X,5F15.6)') (HCOPY(I,J),J=1,NVEC)
    !end do
    !write(u6,*)
    !write(u6,*) ' SCOPY MATRIX:'
    !do I=1,NVEC
    !  write(u6,'(1X,5F15.6)') (SCOPY(I,J),J=1,NVEC)
    !end do
    !write(u6,*)
  end if
  ! --------------------------------------------------------------------
  ! THE UPPER-LEFT NVEC*NVEC SUBMATRICES OF HSMALL AND SSMALL NOW
  ! CONTAINS THE CURRENT HAMILTONIAN AND OVERLAP MATRICES, IN THE
  ! BASIS OF PRESENTLY RETAINED PSI VECTORS. DIAGONALIZE, BUT USE
  ! THE REORDERED MATRICES IN SCOPY, HCOPY. THERE THE BASIS
  ! FUNCTIONS ARE ORDERED BY AGE.
  THR = 1.0e-6_wp
  call SECULAR(MXVEC,NVEC,NRON,HCOPY,SCOPY,VSMALL,ESMALL,SCR,THR)
  ! REORDER THE ELEMENTS OF VSMALL TO GET EIGENVECTORS OF HSMALL. NOTE:
  ! THIS IS NOT THE SAME AS IF WE DIAGONALIZED HSMALL DIRECTLY.
  ! THE DIFFERENCE OCCURS WHENEVER VECTORS ARE THROWN OUT OF THE
  ! CALCULATION IN SECULAR BECAUSE OF LINEAR DEPENDENCE. THE RESULT
  ! WILL DEPEND SLIGHTLY BUT CRITICALLY ON THE ORDER BY WHICH THE
  ! VECTORS WERE ORTHONORMALIZED.
  do I=1,NRON
    do K=1,NVEC
      KK = 1+mod(NVTOT-K,MXVEC)
      SCR(KK) = VSMALL(K,I)
    end do
    do K=1,NVEC
      VSMALL(K,I) = SCR(K)
    end do
  end do
  if (NRON < NRROOT) then
    write(u6,*) 'DIAGRO Error: Linear dependence has reduced'
    write(u6,*) ' the number of solutions to NRON, but you'
    write(u6,*) ' wanted NRROOT soultions.'
    write(u6,'(1X,A,I3)') '  NRON=',NRON
    write(u6,'(1X,A,I3)') 'NRROOT=',NRROOT
    call QUIT(_RC_INTERNAL_ERROR_)
  end if
  ! ORDER THE EIGENFUNCTIONS BY DECREASING OVERLAP WITH THE SPACE
  ! SPANNED BY THE ORIGINALLY SELECTED REFCI ROOTS.
  call DGEMM_('N','N',NVEC,NRON,NVEC,One,PSMALL,MXVEC,VSMALL,MXVEC,Zero,SCR,NVEC)
  II = 1
  do I=1,NRON
    PSEL(I) = DDOT_(NVEC,VSMALL(1,I),1,SCR(II),1)
    II = II+NVEC
  end do
  ! PSEL(I) NOW CONTAINS EXPECTATION VALUE OF PMAT FOR I-TH EIGENVECTOR.
  !write(u6,*) ' ARRAY OF SELECTION AMPLITUDES IN SCR:'
  !write(u6,'(1X,5F15.6)') (PSEL(I),I=1,NRON)
  do I=1,NRON-1
    IMAX = I
    PMAX = PSEL(I)
    do J=I+1,NRON
      if (PSEL(J) >= PMAX) then
        PMAX = PSEL(J)
        IMAX = J
      end if
    end do
    if (IMAX == I) cycle
    PSEL(IMAX) = PSEL(I)
    PSEL(I) = PMAX
    TMP = ESMALL(IMAX)
    ESMALL(IMAX) = ESMALL(I)
    ESMALL(I) = TMP
    do K=1,NVEC
      TMP = VSMALL(K,IMAX)
      VSMALL(K,IMAX) = VSMALL(K,I)
      VSMALL(K,I) = TMP
    end do
  end do
  ! FINALLY, REORDER THE SELECTED ROOTS BY ENERGY:
  do I=1,NRROOT-1
    IMIN = I
    EMIN = ESMALL(I)
    do J=I+1,NRROOT
      if (ESMALL(J) < EMIN) then
        EMIN = ESMALL(J)
        IMIN = J
      end if
    end do
    if (IMIN == I) cycle
    ESMALL(IMIN) = ESMALL(I)
    ESMALL(I) = EMIN
    TMP = PSEL(IMIN)
    PSEL(IMIN) = PSEL(I)
    PSEL(I) = TMP
    do K=1,NVEC
      TMP = VSMALL(K,IMIN)
      VSMALL(K,IMIN) = VSMALL(K,I)
      VSMALL(K,I) = TMP
    end do
  end do
  !write(u6,*) ' EIGENVALUES OF HSMALL. NRON=',NRON
  !write(u6,'(1X,5F15.6)') (ESMALL(I),I=1,NRON)
  !write(u6,*) ' SELECTION WEIGHTS:'
  !write(u6,'(1X,5F15.6)') (   PSEL(I),I=1,NRON)
  !write(u6,*) ' SELECTED EIGENVECTORS:'
  !do I=1,NRROOT
  !  write(u6,'(1X,5F15.6)') (VSMALL(K,I),K=1,NVEC)
  !end do
  !write(u6,*)
  ! --------------------------------------------------------------------
  ! CALCULATE RESIDUAL ARRAYS FOR THE NRROOTS EIGENFUNCTIONS OF HSMALL.
  ! ALSO, USE THE OPPORTUNITY TO FORM MANY OTHER SMALL ARRAYS.
  call mma_allocate(ARR,NRROOT,NRROOT,11,label='ARR')
  call HZLP1(CBUF,SBUF,DBUF,ARR,CSECT,RSECT,XI1,XI2,ICI)
  ! USE THESE SMALLER ARRAYS TO FORM HZERO AND SZERO. THIS IS
  ! OVERLAP AND HAMILTONIAN IN THE BASIS (PSI,RHO,XI1,XI2), WHERE
  ! PSI ARE THE EIGENFUNCTIONS OF HSMALL, RHO ARE RESIDUALS, ETC.
  call HZ(ARR)
  call mma_deallocate(ARR)
  NZ = 4*NRROOT
  !write(u6,*)
  !write(u6,*) ' AFTER HZ CALL. HZERO HAMILTONIAN:'
  !do I=1,NZ
  !  write(u6,'(1X,5F15.6)') (HZERO(I,J),J=1,NZ)
  !end do
  !write(u6,*) ' SZERO:'
  !do I=1,NZ
  !  write(u6,'(1X,5F15.6)') (SZERO(I,J),J=1,NZ)
  !end do
  do I=1,NRROOT
    RNRM(I) = sqrt(SZERO(NRROOT+I,NRROOT+I))
    !EPERT(I) = ESMALL(I)-SZERO(3*NRROOT+I,NRROOT+I)
  end do
  !write(u6,*)
  !write(u6,*) ' PERTURBATION ESTIMATES TO ENERGY:'
  !write(u6,'(1X,5F15.6)') (ESHIFT+EPERT(I),I=1,NRROOT)
  ! --------------------------------------------------------------------
  NCONV = 0
  call TIMING(CPTNOW,DUM1,DUM2,DUM3)
  CPTIT = CPTNOW-CPTOLD
  CPTOLD = CPTNOW
  CPTOT = CPTNOW-CPTSTA
  if (ITER == 1) then
    EDISP = ESMALL(1)+ESHIFT
    write(u6,1234) ITER,NVEC,EDISP,RNRM(1),PSEL(1),CPTIT,CPTOT
  else
    ELOW = ESMALL(1)-ELAST(1)
    if ((ELOW < Zero) .and. (abs(ELOW) <= ETHRE)) NCONV = 1
    EDISP = ESMALL(1)+ESHIFT
    write(u6,1235) ITER,NVEC,EDISP,ELOW,RNRM(1),PSEL(1),CPTIT,CPTOT
  end if
  if (NRROOT > 1) then
    do I=2,NRROOT
      EDISP = ESMALL(I)+ESHIFT
      if (ITER == 1) then
        write(u6,1236) EDISP,RNRM(I),PSEL(I)
      else
        ELOW = ESMALL(I)-ELAST(I)
        if ((ELOW < Zero) .and. (abs(ELOW) <= ETHRE)) NCONV = NCONV+1
        write(u6,1237) EDISP,ELOW,RNRM(I),PSEL(I)
      end if
    end do
    write(u6,*)
  end if
  do I=1,NRROOT
    ELAST(I) = ESMALL(I)
  end do
  if (NCONV == NRROOT) then
    write(u6,*) ' CONVERGENCE IN ENERGY.'
    exit
  end if
  ! --------------------------------------------------------------------
  THR = 1.0e-6_wp
  call SECULAR(MXZ,NZ,NRON,HZERO,SZERO,VZERO,EZERO,SCR,THR)
  !write(u6,*) ' AFTER SECULAR CALL. NRON=',NRON
  !write(u6,*) ' EIGENVALUES & -VECTORS:'
  !do I=1,NRON
  !  write(u6,'(1X,5F15.6)') EZERO(I)
  !  write(u6,'(1X,5F15.6)') (VZERO(K,I),K=1,NZ)
  !end do
  ! ORDER THE EIGENFUNCTIONS BY DECREASING SIZE OF PSI PART.
  call DGEMM_('T','N',NRON,NRROOT,NZ,One,VZERO,MXZ,SZERO,MXZ,Zero,SCR(1+NRON),NRON)
  do I=1,NRON
    II = I
    RSUM = Zero
    do K=1,NRROOT
      II = II+NRON
      RSUM = RSUM+SCR(II)**2
    end do
    SCR(I) = RSUM
  end do
  !write(u6,*)
  !write(u6,*) ' SELECTION CRITERION VECTOR, BEFORE ORDERING:'
  !write(u6,'(1X,5F15.6)') (SCR(I),I=1,NRON)
  do I=1,NRON-1
    IMAX = I
    PMAX = SCR(I)
    do J=I+1,NRON
      if (SCR(J) >= PMAX) then
        PMAX = SCR(J)
        IMAX = J
      end if
    end do
    if (IMAX == I) cycle
    SCR(IMAX) = SCR(I)
    SCR(I) = PMAX
    TMP = EZERO(IMAX)
    EZERO(IMAX) = EZERO(I)
    EZERO(I) = TMP
    do K=1,NZ
      TMP = VZERO(K,IMAX)
      VZERO(K,IMAX) = VZERO(K,I)
      VZERO(K,I) = TMP
    end do
  end do
  !PAM 94-10-30, must reorder as before:
  ! REORDER THE SELECTED ROOTS BY ENERGY:
  do I=1,NRROOT-1
    IMIN = I
    EMIN = EZERO(I)
    do J=I+1,NRROOT
      if (EZERO(J) < EMIN) then
        EMIN = EZERO(J)
        IMIN = J
      end if
    end do
    if (IMIN == I) cycle
    EZERO(IMIN) = EZERO(I)
    EZERO(I) = EMIN
    TMP = SCR(IMIN)
    SCR(IMIN) = SCR(I)
    SCR(I) = TMP
    do K=1,NZ
      TMP = VZERO(K,IMIN)
      VZERO(K,IMIN) = VZERO(K,I)
      VZERO(K,I) = TMP
    end do
  end do
  !PAM 94-10-30, end of update.
  ! NOTE: IF THE UPDATE PART IS SMALL ENOUGH FOR ALL THE FIRST NRROOT
  ! ARRAY, THE CALCULATION HAS CONVERGED.
  NNEW = 0
  !write(u6,*) ' CONVERGENCE CRITERION: SIZE OF UPDATE PART.'
  do I=1,NRROOT
    SQNRM = One-SCR(I)
    !write(u6,*) ' ROOT NR, SQNRM:',I,SQNRM
    if (SQNRM >= SQNLIM) NNEW = NNEW+1
  end do
  !write(u6,*)
  !write(u6,*) ' EIGENVALUES OF THE HZERO HAMILTONIAN:'
  !write(u6,'(1X,5F15.6)') (EZERO(I),I=1,NRON)
  !write(u6,*) ' SELECTION WEIGHTS:'
  !write(u6,'(1X,5F15.6)') (SCR(I),I=1,NRON)
  !write(u6,*) ' EIGENVECTORS:'
  !do I=1,NRON
  !  write(u6,'(1X,5F15.6)') (VZERO(K,I),K=1,NZ)
  !end do
  !write(u6,*) ' NR OF NEW VECTORS SELECTED, NNEW:',NNEW
  if (NNEW == 0) then
    write(u6,*) ' CONVERGENCE IN NORM.'
    exit
  end if
  ! NOTE: A CHANGE HERE. ALWAYS USE ALL THE NRROOT UPDATED VECTORS TO
  ! AVOID OVERWRITING AN EARLY CONVERGED VECTOR (WHICH HAS NEVER BEEN
  ! OUTDATED BY A LATER) BY A VECTOR BELONGING TO ANOTHER ROOT.
  NNEW = NRROOT
  ! --------------------------------------------------------------------
  ! FORM NEW UPDATED VECTORS: SKIP THE FIRST NRROOT-NNEW VECTORS,
  ! WHICH MAKE NO ESSENTIAL IMPROVEMENT.
  !write(u6,*) ' RESET VZERO TO (0,0,0,1) FOR CONVENTIONAL DAVIDSON.'
  !VZERO(:,1:NRROOT) = Zero
  !call DCOPY_(NRROOT,[One],0,VZERO(3*NRROOT+1,1),MXZ+1)
  call HZLP2(CBUF,SBUF,DBUF,CSECT,RSECT,XI1,XI2,CNEW,ICI)
  if (ITER >= MAXIT) then
    write(u6,*) ' UNCONVERGED.'
    exit
  end if
end do
call mma_deallocate(IDC)
call mma_deallocate(IDS)
call mma_deallocate(ELAST)
call mma_deallocate(HCOPY)
call mma_deallocate(PCOPY)
call mma_deallocate(SCOPY)
call mma_deallocate(PSEL)
call mma_deallocate(RNRM)
call mma_deallocate(HSMALL)
call mma_deallocate(SSMALL)
call mma_deallocate(PSMALL)
call mma_deallocate(EZERO)
call mma_deallocate(HZERO)
call mma_deallocate(SZERO)
call mma_deallocate(VZERO)
write(u6,*) ' ',repeat('*',70)
! WRITE CI VECTORS TO LUREST -- CI RESTART FILE.
IDREST = 0
do I=1,NRROOT
  IVEC = 1+mod(NVTOT-NRROOT+I-1,MXVEC)
  IDISK = IDISKC(IVEC)
  do ISTA=1,NCONF,MBUF
    NN = min(MBUF,NCONF+1-ISTA)
    call iDAFILE(LUEIG,2,ICI,NN,IDISK)
    call UPKVEC(NN,ICI,CI(ISTA))
  end do
  call CSFTRA(' CSF',CI,AREF)
  C2REF = Zero
  do IR=1,NREF
    ICSF = IREFX(IR)
    C = CI(ICSF)
    C2REF = C2REF+C**2
  end do
  IR = IROOT(I)
  ECI = ESMALL(I)+ESHIFT
  ENREF = ECI-EREF(IR)
  C2NREF = One-C2REF
  ! WRITE ENERGIES TO PRINTED OUTPUT, AND SAVE TOTAL ENERGIES TO ENGY
  ! FOR LATER PRINTOUT WITH PROPERTIES:
  write(u6,'(A,I3)') '               FINAL RESULTS FOR STATE NR ',I
  write(u6,'(A,I3)') ' CORRESPONDING ROOT OF REFERENCE CI IS NR:',IR
  write(u6,'(A,F15.8)') '            REFERENCE CI ENERGY:',EREF(IR)
  write(u6,'(A,F15.8)') '         EXTRA-REFERENCE WEIGHT:',C2NREF
  if (ICPF == 1) then
    write(u6,'(A,F15.8)') '        ACPF CORRELATION ENERGY:',ENREF
    write(u6,'(A,F15.8)') '                    ACPF ENERGY:',ECI
    ENGY(I,1) = ECI
    ENGY(I,2) = Zero
    ENGY(I,3) = Zero
    call Add_Info('E_MRACPF',[ECI],1,8)
  else
    write(u6,'(A,F15.8)') '          CI CORRELATION ENERGY:',ENREF
    write(u6,'(A,F15.8)') '                      CI ENERGY:',ECI
    ! APPROXIMATE CORRECTIONS FOR UNLINKED QUADRUPLES:
    QDAV = ENREF*C2NREF/C2REF
    EDAV = ECI+QDAV
    QACPF = ENREF*(C2NREF*(One-GFAC))/(C2REF+GFAC*C2NREF)
    EACPF = ECI+QACPF
    write(u6,'(A,F15.8)') '            DAVIDSON CORRECTION:',QDAV
    write(u6,'(A,F15.8)') '               CORRECTED ENERGY:',EDAV
    write(u6,'(A,F15.8)') '                ACPF CORRECTION:',QACPF
    write(u6,'(A,F15.8)') '               CORRECTED ENERGY:',EACPF
    ENGY(I,1) = ECI
    ENGY(I,2) = QDAV
    ENGY(I,3) = QACPF
    call Add_Info('E_MRSDCI',[ECI],1,8)
  end if
  write(u6,*)
  call PRWF_MRCI(CSPCK,INTSY,INDX,CI,JREFX)
  write(u6,*) ' ',repeat('*',70)
  call dDAFILE(LUREST,1,CI,NCONF,IDREST)
end do

call mma_deallocate(CBUF)
call mma_deallocate(SBUF)
call mma_deallocate(DBUF)
call mma_deallocate(CSECT)
call mma_deallocate(RSECT)
call mma_deallocate(XI1)
call mma_deallocate(XI2)
call mma_deallocate(CNEW)
call mma_deallocate(scr)

return

1234 format(1X,I4,1X,I4,1X,F15.8,9X,ES9.2,1X,F6.3,2(1X,F7.1))
1235 format(1X,I4,1X,I4,1X,F15.8,ES9.2,ES9.2,1X,F6.3,2(1X,F7.1))
1236 format(11X,F15.8,9X,ES9.2,1X,F6.3)
1237 format(11X,F15.8,ES9.2,ES9.2,1X,F6.3)

end subroutine MQCT
