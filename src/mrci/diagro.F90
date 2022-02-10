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

subroutine DIAGRO(CI,SGM,CBUF,SBUF,DBUF,AREF,EREF,CSECT,RSECT,XI1,XI2,CNEW,SCR,ICI)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "warnings.h"
#include "mrci.fh"
#include "WrkSpc.fh"
dimension CI(NCONF), SGM(NCONF)
dimension CBUF(MBUF,MXVEC), SBUF(MBUF,MXVEC), DBUF(MBUF)
dimension AREF(NREF,NREF), EREF(NREF)
dimension CSECT(NSECT,MXVEC), RSECT(NSECT,MXVEC)
dimension XI1(NSECT,NRROOT), XI2(NSECT,NRROOT)
dimension CNEW(NSECT,NRROOT), SCR(*)
dimension ICI(MBUF)
dimension IDC(MXVEC), IDS(MXVEC)
dimension HCOPY(MXVEC,MXVEC), SCOPY(MXVEC,MXVEC)
dimension PCOPY(MXVEC,MXVEC)
dimension ELAST(MXROOT), PSEL(MXVEC), RNRM(MXROOT)
!dimension EPERT(MXROOT)

write(6,*)
write(6,*) ('-',I=1,60)
if (ICPF == 0) then
  write(6,*) '   MR SDCI CALCULATION.'
else
  write(6,*) '   MR ACPF CALCULATION.'
end if
write(6,*) ('-',I=1,60)
write(6,*)
write(6,*) '         CONVERGENCE STATISTICS:'
write(6,'(1X,A)') 'ITER NVEC     ENERGIES    LOWERING RESIDUAL SEL.WGT CPU(S) CPU TOT'
ITER = 0
call SETTIM()
call TIMING(CPTNOW,DUM,DUM,DUM)
CPTOLD = CPTNOW
CPTSTA = CPTNOW
! LOOP HEAD FOR CI ITERATIONS:
1000 continue
ITER = ITER+1
! ----------------------------------------------------------------------
! CALCULATE SIGMA ARRAYS FOR SHIFTED HAMILTONIAN IN MCSF BASIS:
do I=1,NNEW
  IVEC = 1+mod(NVTOT-NNEW+I-1,MXVEC)
  IDISK = IDISKC(IVEC)
  do ISTA=1,NCONF,MBUF
    NN = min(MBUF,NCONF+1-ISTA)
    call iDAFILE(LUEIG,2,ICI,NN,IDISK)
    call UPKVEC(NN,ICI,CI(ISTA))
  end do
  call GETMEM('BMN','ALLO','REAL',LBMN,NBMN)
  call GETMEM('IBMN','ALLO','INTE',LIBMN,NBMN)
  call GETMEM('BIAC2','ALLO','REAL',LBIAC2,ISMAX)
  call GETMEM('BICA2','ALLO','REAL',LBICA2,ISMAX)
  call GETMEM('BFIN3','ALLO','REAL',LBFIN3,KBUFF1)
  call GETMEM('AC1','ALLO','REAL',LAC1,ISMAX)
  call GETMEM('AC2','ALLO','REAL',LAC2,ISMAX)
  call GETMEM('BFIN4','ALLO','REAL',LBFIN4,KBUFF1)
  call GETMEM('ABIJ','ALLO','REAL',LABIJ,NVSQ)
  call GETMEM('AIBJ','ALLO','REAL',LAIBJ,NVSQ)
  call GETMEM('AJBI','ALLO','REAL',LAJBI,NVSQ)
  call GETMEM('ASCR1','ALLO','REAL',LASCR1,NVMAX**2)
  call GETMEM('BSCR1','ALLO','REAL',LBSCR1,NVMAX**2)
  call GETMEM('FSCR1','ALLO','REAL',LFSCR1,NVSQ)
  call GETMEM('FSEC','ALLO','REAL',LFSEC,2*NVSQ)
  call GETMEM('BFIN5','ALLO','REAL',LBFIN5,KBUFF1)
  call GETMEM('ASCR2','ALLO','REAL',LASCR2,NVMAX**2)
  call GETMEM('BSCR2','ALLO','REAL',LBSCR2,NVMAX**2)
  call GETMEM('FSCR2','ALLO','REAL',LFSCR2,NVSQ)
  call GETMEM('DBK','ALLO','REAL',LDBK,2*NVSQ)
  !vv call SIGMA(HWORK,CI,SGM)
  !pam call SIGMA(HWORK)
  !PAM04 call SIGMA(HWork(LSGM),HWork(LAREF),HWork(LCI),HWork(LINTSY),HWork(LINDX),HWork(LBMN),HWork(LIBMN),HWork(LBIAC2), &
  !PAM04            HWork(LBICA2),HWork(LBFIN3),HWork(LFIJKL),HWork(LISAB),HWork(LAC1),HWork(LAC2),HWork(LBFIN4),HWork(LABIJ), &
  !PAM04            HWork(LAIBJ),HWork(LAJBI),HWork(LBFIN1),HWork(LASCR1),HWork(LBSCR1),HWork(LFSCR1),HWork(LFSEC),HWork(LFOCK), &
  !PAM04            HWork(LFSCR2),HWork(LDBK),HWork(LCSPCK))
  call SIGMA(SGM,AREF,CI,IWork(LINTSY),IWork(LINDX),Work(LBMN),IWork(LIBMN),Work(LBIAC2),Work(LBICA2),Work(LBFIN3),IWork(LISAB), &
             Work(LAC1),Work(LAC2),Work(LBFIN4),Work(LABIJ),Work(LAIBJ),Work(LAJBI),Work(LASCR1),Work(LBSCR1),Work(LFSCR1), &
             Work(LFSEC),Work(LFOCK),Work(LBFIN5),Work(LASCR2),Work(LBSCR2),Work(LFSCR2),Work(LDBK),IWork(LCSPCK))
  call GETMEM('BFIN5','FREE','REAL',LBFIN5,KBUFF1)
  call GETMEM('ASCR2','FREE','REAL',LASCR2,NVMAX**2)
  call GETMEM('BSCR2','FREE','REAL',LBSCR2,NVMAX**2)
  call GETMEM('FSCR2','FREE','REAL',LFSCR2,NVSQ)
  call GETMEM('DBK','FREE','REAL',LDBK,2*NVSQ)
  call GETMEM('ABIJ','FREE','REAL',LABIJ,NVSQ)
  call GETMEM('AIBJ','FREE','REAL',LAIBJ,NVSQ)
  call GETMEM('AJBI','FREE','REAL',LAJBI,NVSQ)
  call GETMEM('ASCR1','FREE','REAL',LASCR1,NVMAX**2)
  call GETMEM('BSCR1','FREE','REAL',LBSCR1,NVMAX**2)
  call GETMEM('FSCR1','FREE','REAL',LFSCR1,NVSQ)
  call GETMEM('FSEC','FREE','REAL',LFSEC,2*NVSQ)
  call GETMEM('BMN','FREE','REAL',LBMN,NBMN)
  call GETMEM('IBMN','FREE','INTE',LIBMN,NBMN)
  call GETMEM('BIAC2','FREE','REAL',LBIAC2,ISMAX)
  call GETMEM('BICA2','FREE','REAL',LBICA2,ISMAX)
  call GETMEM('BFIN3','FREE','REAL',LBFIN3,KBUFF1)
  call GETMEM('AC1','FREE','REAL',LAC1,ISMAX)
  call GETMEM('AC2','FREE','REAL',LAC2,ISMAX)
  call GETMEM('BFIN4','FREE','REAL',LBFIN4,KBUFF1)
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
! ----------------------------------------------------------------------
! NR OF VECTORS PRESENTLY RETAINED:
NVEC = min(MXVEC,NVTOT)
! NR OF OLD VECTORS RETAINED:
NOLD = NVEC-NNEW
! ----------------------------------------------------------------------
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
do K=1,NVEC
  do L=1,NNEW
    HCOPY(K,L) = 0.0d00
    SCOPY(K,L) = 0.0d00
    PCOPY(K,L) = 0.0d00
  end do
end do
! THEN LOOP OVER BUFFERS. FIRST GET COPIES OF DISK ADDRESSES:
do K=1,NVEC
  IDC(K) = IDISKC(K)
  IDS(K) = IDISKS(K)
end do
do ISTA=1,NCONF,MBUF
  IEND = min(NCONF,ISTA+MBUF-1)
  IBUF = 1+IEND-ISTA
  do K=1,NVEC
    KK = 1+mod(NVTOT-K,MXVEC)
    call iDAFILE(LUEIG,2,ICI,IBUF,IDC(KK))
    call UPKVEC(IBUF,ICI,CBUF(1,K))
    if (K > NNEW) goto 220
    call dDAFILE(LUEIG,2,SBUF(1,K),IBUF,IDS(KK))
220 continue
  end do
  ! --------------------------------------------------------------------
  ! NOTE: AT THIS POINT, THE COLUMNS NR 1..NVEC OF CBUF WILL
  ! CONTAIN THE BUFFERS OF, FIRST, THE NNEW NEWEST PSI ARRAYS,
  ! THEN, THE NOLD ONES FROM EARLIER ITERATIONS.
  ! THE COLUMNS 1..NNEW OF SBUF WILL CONTAIN THE NEWEST NNEW
  ! SIGMA ARRAYS. LEADING DIMENSION OF CBUF AND SBUF IS MBUF. ACTUAL
  ! BUFFER SIZE IS IBUF, WHICH CAN BE SMALLER. ACCUMULATE:
  call DGEMM_('T','N',NVEC,NNEW,IBUF,1.0d0,CBUF,MBUF,CBUF,MBUF,0.0d0,SCR,NVEC)
  KL = 0
  do L=1,NNEW
    do K=1,NVEC
      KL = KL+1
      SCOPY(K,L) = SCOPY(K,L)+SCR(KL)
    end do
  end do
  call DGEMM_('T','N',NVEC,NNEW,IBUF,1.0d0,CBUF,MBUF,SBUF,MBUF,0.0d0,SCR,NVEC)
  KL = 0
  do L=1,NNEW
    do K=1,NVEC
      KL = KL+1
      HCOPY(K,L) = HCOPY(K,L)+SCR(KL)
    end do
  end do
! ALSO, UPDATE PSMALL, WHICH IS USED FOR SELECTION.
  if (ISTA > IREFX(NRROOT)) goto 200
  do I=1,NRROOT
    IR = IROOT(I)
    IRR = IREFX(IR)
    if (IRR < ISTA) goto 250
    if (IRR > IEND) goto 250
    IPOS = IRR+1-ISTA
    do L=1,NNEW
      do K=1,NVEC
        PCOPY(K,L) = PCOPY(K,L)+CBUF(IPOS,K)*CBUF(IPOS,L)
      end do
    end do
250 continue
  end do
200 continue
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
  write(6,*)
  write(6,*) ' HSMALL MATRIX:'
  do I=1,NVEC
    write(6,'(1X,5F15.6)') (HSMALL(I,J),J=1,NVEC)
  end do
  write(6,*)
  write(6,*) ' SSMALL MATRIX:'
  do I=1,NVEC
    write(6,'(1X,5F15.6)') (SSMALL(I,J),J=1,NVEC)
  end do
  write(6,*)
  write(6,*) ' PSMALL MATRIX:'
  do I=1,NVEC
    write(6,'(1X,5F15.6)') (PSMALL(I,J),J=1,NVEC)
  end do
  !write(6,*)
  !write(6,*)
  !write(6,*) ' HCOPY MATRIX:'
  !do I=1,NVEC
  !  write(6,'(1X,5F15.6)') (HCOPY(I,J),J=1,NVEC)
  !end do
  !write(6,*)
  !write(6,*) ' SCOPY MATRIX:'
  !do I=1,NVEC
  !  write(6,'(1X,5F15.6)') (SCOPY(I,J),J=1,NVEC)
  !end do
  !write(6,*)
end if
! -------------------------------------------------------------------
! THE UPPER-LEFT NVEC*NVEC SUBMATRICES OF HSMALL AND SSMALL NOW
! CONTAINS THE CURRENT HAMILTONIAN AND OVERLAP MATRICES, IN THE
! BASIS OF PRESENTLY RETAINED PSI VECTORS. DIAGONALIZE, BUT USE
! THE REORDERED MATRICES IN SCOPY, HCOPY,DCOPY. THERE THE BASIS
! FUNCTIONS ARE ORDERED BY AGE.
THR = 1.0D-06
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
  write(6,*) 'DIAGRO Error: Linear dependence has reduced'
  write(6,*) ' the number of solutions to NRON, but you'
  write(6,*) ' wanted NRROOT soultions.'
  write(6,'(1X,A,I3)') '  NRON=',NRON
  write(6,'(1X,A,I3)') 'NRROOT=',NRROOT
  call QUIT(_RC_INTERNAL_ERROR_)
end if
! ORDER THE EIGENFUNCTIONS BY DECREASING OVERLAP WITH THE SPACE
! SPANNED BY THE ORIGINALLY SELECTED REFCI ROOTS.
call DGEMM_('N','N',NVEC,NRON,NVEC,1.0d0,PSMALL,MXVEC,VSMALL,MXVEC,0.0d0,SCR,NVEC)
II = 1
do I=1,NRON
  PSEL(I) = DDOT_(NVEC,VSMALL(1,I),1,SCR(II),1)
  II = II+NVEC
end do
! PSEL(I) NOW CONTAINS EXPECTATION VALUE OF PMAT FOR I-TH EIGENVECTOR.
!write(6,*) ' ARRAY OF SELECTION AMPLITUDES IN SCR:'
!write(6,'(1X,5F15.6)') (PSEL(I),I=1,NRON)
do I=1,NRON-1
  IMAX = I
  PMAX = PSEL(I)
  do J=I+1,NRON
    if (PSEL(J) < PMAX) goto 360
    PMAX = PSEL(J)
    IMAX = J
360 continue
  end do
  if (IMAX == I) goto 380
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
380 continue
end do
! FINALLY, REORDER THE SELECTED ROOTS BY ENERGY:
do I=1,NRROOT-1
  IMIN = I
  EMIN = ESMALL(I)
  do J=I+1,NRROOT
    if (ESMALL(J) >= EMIN) goto 1360
    EMIN = ESMALL(J)
    IMIN = J
1360 continue
  end do
  if (IMIN == I) goto 1380
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
1380 continue
end do
!write(6,*) ' EIGENVALUES OF HSMALL. NRON=',NRON
!write(6,'(1X,5F15.6)') (ESMALL(I),I=1,NRON)
!write(6,*) ' SELECTION WEIGHTS:'
!write(6,'(1X,5F15.6)') (   PSEL(I),I=1,NRON)
!write(6,*) ' SELECTED EIGENVECTORS:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (VSMALL(K,I),K=1,NVEC)
!end do
!write(6,*)
! ----------------------------------------------------------------------
! CALCULATE RESIDUAL ARRAYS FOR THE NRROOTS EIGENFUNCTIONS OF HSMALL.
! ALSO, USE THE OPPORTUNITY TO FORM MANY OTHER SMALL ARRAYS.
call GETMEM('ARR','ALLO','REAL',LARR,11*NRROOT**2)
call HZLP1(CBUF,SBUF,DBUF,WORK(LARR),CSECT,RSECT,XI1,XI2,ICI)
! USE THESE SMALLER ARRAYS TO FORM HZERO AND SZERO. THIS IS
! OVERLAP AND HAMILTONIAN IN THE BASIS (PSI,RHO,XI1,XI2), WHERE
! PSI ARE THE EIGENFUNCTIONS OF HSMALL, RHO ARE RESIDUALS, ETC.
call HZ(WORK(LARR))
call GETMEM('ARR','FREE','REAL',LARR,11*NRROOT**2)
NZ = 4*NRROOT
!write(6,*)
!write(6,*) ' AFTER HZ CALL. HZERO HAMILTONIAN:'
!do I=1,NZ
!  write(6,'(1X,5F15.6)') (HZERO(I,J),J=1,NZ)
!end do
!write(6,*) ' SZERO:'
!do I=1,NZ
!  write(6,'(1X,5F15.6)') (SZERO(I,J),J=1,NZ)
!end do
do I=1,NRROOT
  RNRM(I) = sqrt(SZERO(NRROOT+I,NRROOT+I))
  !EPERT(I) = ESMALL(I)-SZERO(3*NRROOT+I,NRROOT+I)
end do
!write(6,*)
!write(6,*) ' PERTURBATION ESTIMATES TO ENERGY:'
!write(6,'(1X,5F15.6)') (ESHIFT+EPERT(I),I=1,NRROOT)
! ----------------------------------------------------------------------
NCONV = 0
call TIMING(CPTNOW,DUM,DUM,DUM)
CPTIT = CPTNOW-CPTOLD
CPTOLD = CPTNOW
CPTOT = CPTNOW-CPTSTA
if (ITER == 1) then
  EDISP = ESMALL(1)+ESHIFT
  write(6,1234) ITER,NVEC,EDISP,RNRM(1),PSEL(1),CPTIT,CPTOT
else
  ELOW = ESMALL(1)-ELAST(1)
  if ((ELOW < 0.0d00) .and. (abs(ELOW) <= ETHRE)) NCONV = 1
  EDISP = ESMALL(1)+ESHIFT
  write(6,1235) ITER,NVEC,EDISP,ELOW,RNRM(1),PSEL(1),CPTIT,CPTOT
end if
if (NRROOT > 1) then
  do I=2,NRROOT
    EDISP = ESMALL(I)+ESHIFT
    if (ITER == 1) then
      write(6,1236) EDISP,RNRM(I),PSEL(I)
    else
      ELOW = ESMALL(I)-ELAST(I)
      if ((ELOW < 0.0d00) .and. (abs(ELOW) <= ETHRE)) NCONV = NCONV+1
      write(6,1237) EDISP,ELOW,RNRM(I),PSEL(I)
    end if
  end do
  write(6,*)
end if
do I=1,NRROOT
  ELAST(I) = ESMALL(I)
end do
1234 format(1X,I4,1X,I4,1X,F15.8,9X,D9.2,1X,F6.3,2(1X,F7.1))
1235 format(1X,I4,1X,I4,1X,F15.8,D9.2,D9.2,1X,F6.3,2(1X,F7.1))
1236 format(11X,F15.8,9X,D9.2,1X,F6.3)
1237 format(11X,F15.8,D9.2,D9.2,1X,F6.3)
if (NCONV == NRROOT) then
  write(6,*) ' CONVERGENCE IN ENERGY.'
  goto 2000
end if
! ------------------------------------------------------------------
THR = 1.0D-06
call SECULAR(MXZ,NZ,NRON,HZERO,SZERO,VZERO,EZERO,SCR,THR)
!write(6,*) ' AFTER SECULAR CALL. NRON=',NRON
!write(6,*) ' EIGENVALUES & -VECTORS:'
!do I=1,NRON
!  write(6,'(1X,5F15.6)') EZERO(I)
!  write(6,'(1X,5F15.6)') (VZERO(K,I),K=1,NZ)
!end do
! ORDER THE EIGENFUNCTIONS BY DECREASING SIZE OF PSI PART.
call DGEMM_('T','N',NRON,NRROOT,NZ,1.0d0,VZERO,MXZ,SZERO,MXZ,0.0d0,SCR(1+NRON),NRON)
do I=1,NRON
  II = I
  SUM = 0.0d00
  do K=1,NRROOT
    II = II+NRON
    SUM = SUM+SCR(II)**2
  end do
  SCR(I) = SUM
end do
!write(6,*)
!write(6,*) ' SELECTION CRITERION VECTOR, BEFORE ORDERING:'
!write(6,'(1X,5F15.6)') (SCR(I),I=1,NRON)
do I=1,NRON-1
  IMAX = I
  PMAX = SCR(I)
  do J=I+1,NRON
    if (SCR(J) < PMAX) goto 460
    PMAX = SCR(J)
    IMAX = J
460 continue
  end do
  if (IMAX == I) goto 480
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
480 continue
end do
!PAM 94-10-30, must reorder as before:
! REORDER THE SELECTED ROOTS BY ENERGY:
do I=1,NRROOT-1
  IMIN = I
  EMIN = EZERO(I)
  do J=I+1,NRROOT
    if (EZERO(J) >= EMIN) goto 1460
    EMIN = EZERO(J)
    IMIN = J
1460 continue
  end do
  if (IMIN == I) goto 1480
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
1480 continue
end do
!PAM 94-10-30, end of update.
! NOTE: IF THE UPDATE PART IS SMALL ENOUGH FOR ALL THE FIRST NRROOT
! ARRAY, THE CALCULATION HAS CONVERGED.
NNEW = 0
!write(6,*) ' CONVERGENCE CRITERION: SIZE OF UPDATE PART.'
do I=1,NRROOT
  SQNRM = 1.0d00-SCR(I)
  !write(6,*) ' ROOT NR, SQNRM:',I,SQNRM
  if (SQNRM < SQNLIM) goto 490
  NNEW = NNEW+1
490 continue
end do
!write(6,*)
!write(6,*) ' EIGENVALUES OF THE HZERO HAMILTONIAN:'
!write(6,'(1X,5F15.6)') (EZERO(I),I=1,NRON)
!write(6,*) ' SELECTION WEIGHTS:'
!write(6,'(1X,5F15.6)') (SCR(I),I=1,NRON)
!write(6,*) ' EIGENVECTORS:'
!do I=1,NRON
!  write(6,'(1X,5F15.6)') (VZERO(K,I),K=1,NZ)
!end do
!write(6,*) ' NR OF NEW VECTORS SELECTED, NNEW:',NNEW
if (NNEW == 0) then
  write(6,*) ' CONVERGENCE IN NORM.'
  goto 2000
end if
! NOTE: A CHANGE HERE. ALWAYS USE ALL THE NRROOT UPDATED VECTORS TO
! AVOID OVERWRITING AN EARLY CONVERGED VECTOR (WHICH HAS NEVER BEEN
! OUTDATED BY A LATER) BY A VECTOR BELONGING TO ANOTHER ROOT.
NNEW = NRROOT
! ----------------------------------------------------------------------
! FORM NEW UPDATED VECTORS: SKIP THE FIRST NRROOT-NNEW VECTORS,
! WHICH MAKE NO ESSENTIAL IMPROVEMENT.
!write(6,*) ' RESET VZERO TO (0,0,0,1) FOR CONVENTIONAL DAVIDSON.'
!call DCOPY_(NRROOT*MXZ,[0.0D00],0,VZERO,1)
!call DCOPY_(NRROOT,[1.0D00],0,VZERO(3*NRROOT+1,1),MXZ+1)
call HZLP2(CBUF,SBUF,DBUF,CSECT,RSECT,XI1,XI2,CNEW,ICI)
if (ITER < MAXIT) goto 1000
write(6,*) ' UNCONVERGED.'
2000 continue
write(6,*) ' ',('*',III=1,70)
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
  C2REF = 0.0d00
  do IR=1,NREF
    ICSF = IREFX(IR)
    C = CI(ICSF)
    C2REF = C2REF+C**2
  end do
  IR = IROOT(I)
  ECI = ESMALL(I)+ESHIFT
  ENREF = ECI-EREF(IR)
  C2NREF = 1.0d00-C2REF
  ! WRITE ENERGIES TO PRINTED OUTPUT, AND SAVE TOTAL ENERGIES TO ENGY
  ! FOR LATER PRINTOUT WITH PROPERTIES:
  write(6,'(A,I3)') '               FINAL RESULTS FOR STATE NR ',I
  write(6,'(A,I3)') ' CORRESPONDING ROOT OF REFERENCE CI IS NR:',IR
  write(6,'(A,F15.8)') '            REFERENCE CI ENERGY:',EREF(IR)
  write(6,'(A,F15.8)') '         EXTRA-REFERENCE WEIGHT:',C2NREF
  if (ICPF == 1) then
    write(6,'(A,F15.8)') '        ACPF CORRELATION ENERGY:',ENREF
    write(6,'(A,F15.8)') '                    ACPF ENERGY:',ECI
    ENGY(I,1) = ECI
    ENGY(I,2) = 0.0d00
    ENGY(I,3) = 0.0d00
    call Add_Info('E_MRACPF',[ECI],1,8)
  else
    write(6,'(A,F15.8)') '          CI CORRELATION ENERGY:',ENREF
    write(6,'(A,F15.8)') '                      CI ENERGY:',ECI
    ! APPROXIMATE CORRECTIONS FOR UNLINKED QUADRUPLES:
    QDAV = ENREF*C2NREF/C2REF
    EDAV = ECI+QDAV
    QACPF = ENREF*(C2NREF*(1.0d00-GFAC))/(C2REF+GFAC*C2NREF)
    EACPF = ECI+QACPF
    write(6,'(A,F15.8)') '            DAVIDSON CORRECTION:',QDAV
    write(6,'(A,F15.8)') '               CORRECTED ENERGY:',EDAV
    write(6,'(A,F15.8)') '                ACPF CORRECTION:',QACPF
    write(6,'(A,F15.8)') '               CORRECTED ENERGY:',EACPF
    ENGY(I,1) = ECI
    ENGY(I,2) = QDAV
    ENGY(I,3) = QACPF
    call Add_Info('E_MRSDCI',[ECI],1,8)
  end if
  write(6,*)
  !PAM04 call PRWF_MRCI (HWORK(LCSPCK),HWORK(LINTSY),HWORK(LINDX),CI,HWORK(LJREFX) )
  call PRWF_MRCI(IWORK(LCSPCK),IWORK(LINTSY),IWORK(LINDX),CI,IWORK(LJREFX))
  write(6,*) ' ',('*',III=1,70)
  call dDAFILE(LUREST,1,CI,NCONF,IDREST)
end do
call XFlush(6)

return

end subroutine DIAGRO
