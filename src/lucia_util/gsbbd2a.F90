!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Jeppe Olsen                                      *
!               2026, Meng Wang                                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GSBBD2A(RHO2,RHO2S,RHO2A,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,NSB,NCB,MXPNGAS,NOBPTS, &
                   IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,X,NSMOB,SCLFAC,IPACK)
! SUBROUTINE GSBBD2A --> 37
!
! Contributions to two-electron density matrix from column excitations
!
! GAS version, '96, Jeppe Olsen
!
! =====
! Input
! =====
! RHO2  : two body density matrix to be updated
! NACOB : Number of active orbitals
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP : String group of columns
! NROW : Number of rows in S and C block
! NGAS : Number of active spaces
! ISEL : Number of electrons per AS for S block
! ICEL : Number of electrons per AS for C block
! CB   : Input C block
! MXPNGAS : Max number of AS spaces (program parameter)
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB  : Number of symmetries of orbitals
! MAXI   : Largest Number of "spectator strings" treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
! IPACK  : Should we pack the density?
!
! ======
! Output
! ======
! RHO2, RHO2S, RHO2A : Updated density block
!
! =======
! Scratch
! =======
!
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S, I2, XI2S : For holding creations/annihilations type and symmetry
!
! Jeppe Olsen, Fall of 96

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Para_Info, only: MyRank, nProcs
use lucia_runtime, only: LUCIA_OPTIMIZATIONS_ENABLED
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _CUDA_BLAS_
use, intrinsic :: iso_c_binding, only: c_int64_t
use GSBBD2A_CUDA_INTERFACE, only: LUCIA_GSBBD2A_CUDA_BEGIN, LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN, &
                                  LUCIA_GSBBD2A_CUDA_DENSITY_END, LUCIA_GSBBD2A_CUDA_BLOCK_END, &
                                  LUCIA_GSBBD2A_CUDA_END, LUCIA_GSBBD2A_CUDA_ROUTE
#endif
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*)
integer(kind=iwp), intent(in) :: NACOB, ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, NSB, NCB, ISEL(NGAS), ICEL(NGAS), MXPNGAS, &
                                 NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), MAXI, MAXK, NSMOB
real(kind=wp), intent(in) :: SB(NROW,NSB), CB(NROW,NCB), SCLFAC
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(MAXK,*), XI2S(MAXK,*), X(*)
integer(kind=iwp), intent(_OUT_) :: I1(MAXK,*), I2(MAXK,*)
logical(kind=iwp), intent(in) :: IPACK
integer(kind=iwp) :: I, I1IK, I1JL, IBOT, IDXSM, IDXTP, IFIRST, IFRST, II12, IIK, IIKE, IIPART, IJL, IJLE, IKBOFF, IKOBSM, IKOFF, &
                     IKSM, IOFF, ISM, ITOP, ITP(256), ITYP, J, JFRST, JLBOFF, JLOBSM, JLOFF, JLSM, JOFF, JSM, JTP(256), JTYP, K, &
                     K12, KBOT, KEND, KFRST, KOFF, KSM, KTOP, KTP(256), KTYP, L, LDUMMY, LOFF, LSM, LTP(256), LTYP, NDXTP, NI, &
                     NIBTC, NIK, NJ, NJL, NK, NKBTC, NL, NONEW, NPART, NPARTSZ
real(kind=wp) :: FACTOR
#ifdef _CUDA_BLAS_
integer(c_int64_t) :: CudaStatus
logical :: CudaSession, DensitySession, BlockFused
#endif

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ================='
write(u6,*) ' GSBBD2A in action'
write(u6,*) ' ================='
write(u6,*)
write(u6,*) ' Occupation of active left strings'
call IWRTMA(ISEL,1,NGAS,1,NGAS)
write(u6,*) ' Occupation of active Right strings'
call IWRTMA(ICEL,1,NGAS,1,NGAS)
#endif

IFRST = 1
JFRST = 1

!SVC: determine optimum number of partitions as the lowest multiple of
!     NPROCS that satisfies a block size smaller than MAXI:
NPART = 0
do
  NPART = NPART+NPROCS
  NPARTSZ = max(NROW-1,0)/NPART+1
  if (NPARTSZ <= MAXI) exit
end do

! Type of single excitations that connects the two column strings
call DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NGAS,ISEL,ICEL)
! Symmetry of Double excitation that connects IBSM and JBSM
IDXSM = Mul(ISCSM,ICCSM)
if (IDXSM /= 0) then
#ifdef _CUDA_BLAS_
  DensitySession = .false.
  if (LUCIA_OPTIMIZATIONS_ENABLED() .and. NPROCS == 1 .and. NACOB > 0) then
    CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_BEGIN(RHO2,RHO2S,RHO2A,int(NACOB,c_int64_t), &
                                                  merge(1_c_int64_t,0_c_int64_t,IPACK))
    if (CudaStatus == -1_c_int64_t) then
      call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
    else
      DensitySession = CudaStatus == 1_c_int64_t
    end if
  end if
#endif
# ifdef _DEBUGPRINT_
  write(u6,*) ' ISCSM,ICCSM ',ISCSM,ICCSM
# endif
  do IDXTP=1,NDXTP
    ITYP = ITP(IDXTP)
    JTYP = JTP(IDXTP)
    KTYP = KTP(IDXTP)
    LTYP = LTP(IDXTP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ITYP JTYP KTYP LTYP ',ITYP,JTYP,KTYP,LTYP
#   endif
    do IKOBSM=1,NSMOB
      JLOBSM = Mul(IKOBSM,IDXSM)
      if (JLOBSM == 0) cycle
      ! types + symmetries defined => K strings are defined
      KFRST = 1
      ! Loop over of symmetry of i orbitals
      do ISM=1,NSMOB
        KSM = Mul(ISM,IKOBSM)
        NI = NOBPTS(ITYP,ISM)
        NK = NOBPTS(KTYP,KSM)
        if ((NI == 0) .or. (NK == 0)) cycle
        ! Loop over batches of j orbitals
        outer: do JSM=1,NSMOB
          LSM = Mul(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)
          if ((NJ == 0) .or. (NL == 0)) cycle outer

          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)

          if ((IOFF < KOFF) .or. (JOFF < LOFF)) cycle outer

          ! ============================================================
          !                  Use N-2 projection method
          ! ============================================================

          IFIRST = 1
          ! Loop over batches of I strings
          X(1:NI*NJ*NK*NL) = Zero
#ifdef _CUDA_BLAS_
          CudaSession = .false.
          if (LUCIA_OPTIMIZATIONS_ENABLED() .and. NPROCS == 1 .and. NROW > 0 .and. NSB > 0 .and. NCB > 0 .and. NI*NJ*NK*NL > 0) then
            CudaStatus = LUCIA_GSBBD2A_CUDA_BEGIN(X,SB,CB,NROW,NSB,NCB,NI*NJ*NK*NL)
            if (CudaStatus == -1_c_int64_t) then
              call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
            else
              CudaSession = CudaStatus == 1_c_int64_t
              if ((.not. CudaSession) .and. DensitySession) then
                CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_END()
                if (CudaStatus == -1_c_int64_t) then
                  call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
                end if
                DensitySession = .false.
              end if
            end if
          else if (DensitySession) then
            CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_END()
            if (CudaStatus == -1_c_int64_t) then
              call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
            end if
            DensitySession = .false.
          end if
#endif
          do IIPART=1+MYRANK,NPART,NPROCS
            IBOT = 1+(IIPART-1)*NPARTSZ
            ITOP = min(IBOT+NPARTSZ-1,NROW)
            NIBTC = ITOP-IBOT+1
            if (NIBTC <= 0) exit
            ! Loop over batches of intermediate strings
            KBOT = 1-MAXK
            KTOP = 0
            do
              KBOT = KBOT+MAXK
              KTOP = KTOP+MAXK

              ! =======================================================
              !
              ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
              !
              ! =======================================================

              JLBOFF = 1
              if ((JSM == LSM) .and. (JTYP == LTYP)) then
                NJL = nTri_Elem(NJ)
                JLSM = 1
              else
                NJL = NJ*NL
                JLSM = 0
              end if
              ! Obtain all double excitations from this group of K strings
              II12 = 1
              K12 = 1
              call ADADST_GAS(1,JSM,JTYP,NJ,1,LSM,LTYP,NL,ICCTP,ICCSM,IGRP,KBOT,KTOP,I2,XI2S,MAXK,NKBTC,KEND,JFRST,KFRST,II12,K12, &
                              SCLFAC)

              JFRST = 0
              KFRST = 0

              if (NKBTC == 0) then
#ifdef _CUDA_BLAS_
                if (CudaSession) then
                  CudaStatus = LUCIA_GSBBD2A_CUDA_END()
                  if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
                  CudaSession = .false.
                end if
#endif
                cycle outer
              end if

              ! =======================================================
              !
              ! obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
              !
              ! =======================================================

              IKBOFF = 1
              if ((ISM == KSM) .and. (ITYP == KTYP)) then
                NIK = nTri_Elem(NI)
                IKSM = 1
              else
                NIK = NI*NK
                IKSM = 0
              end if
              ! Obtain all double excitations from this group of K strings
              II12 = 2
              K12 = 1
              if (IFRST == 1) KFRST = 1
              call ADADST_GAS(1,ISM,ITYP,NI,1,KSM,KTYP,NK,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,IFRST,KFRST,II12,K12, &
                              One)

              IFRST = 0
              KFRST = 0

              if (NKBTC == 0) then
#ifdef _CUDA_BLAS_
                if (CudaSession) then
                  CudaStatus = LUCIA_GSBBD2A_CUDA_END()
                  if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
                  CudaSession = .false.
                end if
#endif
                cycle outer
              end if

              if (IFIRST == 1) then
                FACTOR = Zero
              else
                FACTOR = One
              end if

#ifdef _CUDA_BLAS_
              CudaStatus = 0_c_int64_t
              if (LUCIA_OPTIMIZATIONS_ENABLED() .and. NPROCS == 1) then
                CudaStatus = LUCIA_GSBBD2A_CUDA_ROUTE(X,SB,CB,I1,XI1S,I2,XI2S,NROW,NSB,NCB,IBOT,NIBTC,MAXK,NKBTC,NI,NK,NJ,NL, &
                                                      NIK,NJL,IKSM,JLSM,FACTOR)
              end if
              if (CudaStatus == -1_c_int64_t) then
                call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
              else if (CudaStatus /= 1_c_int64_t) then
                if (DensitySession) then
                  CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_END()
                  if (CudaStatus == -1_c_int64_t) then
                    call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
                  end if
                  DensitySession = .false.
                end if
#endif
              ! Gather C for the CPU fallback.
              J = 0
              L = 1
              do IJL=1,NJL
                call NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                I1JL = (L-1)*NJ+J
                if (JLSM /= 0) then
                  IJLE = nTri_Elem(J-1)+L
                else
                  IJLE = IJL
                end if
                JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
                if ((JLSM == 1) .and. (J == L)) then
                  CSCR(JLOFF:JLOFF+NKBTC*NIBTC-1) = Zero
                else
                  call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I2(:,I1JL),XI2S(:,I1JL))
                end if
              end do

              ! Gather S for the CPU fallback.
              I = 0
              K = 1
              do IIK=1,NIK
                call NXTIJ(I,K,NI,NK,IKSM,NONEW)
                I1IK = (K-1)*NI+I
                if (IKSM /= 0) then
                  IIKE = nTri_Elem(I-1)+K
                else
                  IIKE = IIK
                end if
                IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
                if ((IKSM == 1) .and. (I == K)) then
                  SSCR(IKOFF:IKOFF+NKBTC*NIBTC-1) = Zero
                else
                  call MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,I1(:,I1IK),XI1S(:,I1IK))
                end if
              end do

              ! ==================================================================
              !
              ! RHO2C(ik,jl)  = RHO2C(ik,jl) - sum(Ia,Kb)SB(Ia,Kb,ik)*CB(Ia,Kb,jl)
              !
              ! ==================================================================

              ! The minus ??
              !
              ! Well, the density matrices are constructed as
              !
              ! <I!a+i a+k aj al!> = -sum(K) <I!a+ia+k!K><J!aj al!K>, and
              ! the latter matrices are the ones we are constructing

              IOFF = IOBPTS(ITYP,ISM)
              JOFF = IOBPTS(JTYP,JSM)
              KOFF = IOBPTS(KTYP,KSM)
              LOFF = IOBPTS(LTYP,LSM)
              LDUMMY = NKBTC*NIBTC
#             ifdef _DEBUGPRINT_
              write(u6,*) ' CSCR matrix'
              call WRTMAT(CSCR,LDUMMY,NJL,LDUMMY,NJL)
              write(u6,*) ' SSCR matrix'
              call WRTMAT(SSCR,LDUMMY,NIK,LDUMMY,NIK)
#             endif

              LDUMMY = NKBTC*NIBTC
              !    MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP)
              call MATML7(X,SSCR,CSCR,NIK,NJL,LDUMMY,NIK,LDUMMY,NJL,FACTOR,-One,1)
#ifdef _CUDA_BLAS_
              end if
#endif
              IFIRST = 0
#             ifdef _DEBUGPRINT_
              write(u6,*) ' Updated X matrix IK,JL,IK,JL',NIK,NJL,NIK,NJL
              call WRTMAT(X,NIK,NJL,NIK,NJL)
#             endif

              if (KEND /= 0) exit
            end do
            ! End of loop over partitionings of resolution strings
          end do
          ! Rho2(ik,jl) has been constructed for ik,jl belonging to
          ! Scatter out to density matrix
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
          !write(u6,*) 'I, J, K, L offsets & IPACK :',IOFF,JOFF,KOFF,LOFF,Ipack
#ifdef _CUDA_BLAS_
          BlockFused = .false.
          if (CudaSession) then
            CudaStatus = LUCIA_GSBBD2A_CUDA_BLOCK_END(int(NI,c_int64_t),int(IOFF,c_int64_t),int(NJ,c_int64_t), &
                                                       int(JOFF,c_int64_t),int(NK,c_int64_t),int(KOFF,c_int64_t), &
                                                       int(NL,c_int64_t),int(LOFF,c_int64_t),int(NACOB,c_int64_t), &
                                                       merge(1_c_int64_t,0_c_int64_t,IPACK))
            if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
            BlockFused = CudaStatus == 1_c_int64_t
            CudaSession = .false.
          else if (DensitySession) then
            CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_END()
            if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
            DensitySession = .false.
          end if
#endif
#ifdef _CUDA_BLAS_
          if (.not. BlockFused) then
#endif
          call ADTOR2(RHO2,RHO2S,RHO2A,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NACOB,IPACK)
#ifdef _CUDA_BLAS_
          end if
#endif
          !    ADTOR2(RHO2,RHO2T,ITYPE,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)

          !write(u6,*) ' updated density matrix A',' norb = 4 ',norb
          !write(u6,*) ' offset ','IOFF,JOFF,KOFF,LOFF',IOFF,JOFF,KOFF,LOFF
          !call prsym(rho2s,nTri_Elem(NORB))

        end do outer
      end do
    end do
  end do
#ifdef _CUDA_BLAS_
  if (DensitySession) then
    CudaStatus = LUCIA_GSBBD2A_CUDA_DENSITY_END()
    if (CudaStatus == -1_c_int64_t) call SYSABENDMSG('lucia_util/gsbbd2a','CUDA execution failed','')
    DensitySession = .false.
  end if
#endif
end if

end subroutine GSBBD2A
