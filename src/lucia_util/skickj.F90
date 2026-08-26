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
! Copyright (C) 1994,1997, Jeppe Olsen                                 *
!               2026, Meng Wang                                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SKICKJ(SKII,CKJJ,NKA,NKB,XIJKL,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,FACS,IROUTE)
! Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
!          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
!
! Jeppe Olsen, Spring of 94
!
! : Note : Route 1 has retired, March 97

use lucia_data, only: MXPTSOB, SKICKJ_TINY_K_MAX, SKICKJ_TINY_N_MAX
#ifdef _CUDA_BLAS_
use, intrinsic :: iso_c_binding, only: c_int64_t
use LUCIA_CUDA_INTERFACE, only: LUCIA_SKICKJ_CUDA_ROUTE3
use Constants, only: One
#endif
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: SKII(*), XIJKL(*)
integer(kind=iwp), intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, KBIB(MAXK,*), KBJB(MAXK,*), IKORD, IROUTE
real(kind=wp), intent(in) :: CKJJ(*), XKBIB(MAXK,*), XKBJB(MAXK,*), FACS
integer(kind=iwp) :: IB, ICOFF, IKINTOF, IMAX, INTOF, ISOFF, JB, JKINTOF, K, KB, KK, L, LL
real(kind=wp) :: FACTOR, SGNK, SGNL, XIJILS(MXPTSOB)
#ifdef _CUDA_BLAS_
integer(kind=c_int64_t) :: CudaStatus, NIBCUDA, NJBCUDA
#endif

! To get rid of annoying and incorrect compiler warnings
JKINTOF = 0
IKINTOF = 0

if ((NI > MXPTSOB) .or. (NJ > MXPTSOB) .or. (NK > MXPTSOB) .or. (NL > MXPTSOB)) then
  write(u6,*) ' SKICKJ : Too many orbs : > MXPTSOB'
  write(u6,*) ' N, MXPTSOB ',max(NI,NJ,NK,NL),MXPTSOB
  !stop ' Redim MXPTSOB'
  call SYSABENDMSG('lucia_util/skickj','Redim MXPTSOB','')
end if

if (IROUTE == 3) then
  ! S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
# ifdef _CUDA_BLAS_
  if ((NKA > 0) .and. (NKB > 0) .and. (NI > 0) .and. (NJ > 0) .and. (NK > 0) .and. (NL > 0) .and. (MAXK > 0) .and. (NKB <= MAXK) &
      .and. ((IKORD == 0) .or. (IKORD == 1)) .and. (FACS == One)) then
    NIBCUDA = maxval(KBIB(1:NKB,1:NK))
    NJBCUDA = maxval(KBJB(1:NKB,1:NL))
    if ((NIBCUDA > 0) .and. (NJBCUDA > 0)) then
      CudaStatus = LUCIA_SKICKJ_CUDA_ROUTE3(SKII,CKJJ,XIJKL,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,NIBCUDA,NJBCUDA, &
                                            FACS)
      if (CudaStatus == 1_c_int64_t) return
    end if
  end if
# endif
  do KB=1,NKB
    ! Number of nonvanishing connections from KB
    LL = 0
    KK = 0
    do L=1,NL
      if (KBJB(KB,L) /= 0) LL = LL+1
    end do
    do K=1,NK
      if (KBIB(KB,K) /= 0) KK = KK+1
    end do

    if ((KK /= 0) .and. (LL /= 0)) then
      do K=1,NK
        IB = KBIB(KB,K)
        if (IB /= 0) then
          SGNK = XKBIB(KB,K)
          do L=1,NL
            JB = KBJB(KB,L)
#           ifdef _DEBUGPRINT_
            write(u6,*) ' KB,K,L,IB,JB',KB,K,L,IB,JB
#           endif
            if (JB /= 0) then
              SGNL = XKBJB(KB,L)
              FACTOR = SGNK*SGNL
              ! We have now a IB and Jb string, let's do it
              ISOFF = (IB-1)*NI*NKA+1
              ICOFF = (JB-1)*NJ*NKA+1
              INTOF = ((L-1)*NK+K-1)*NI*NJ+1
              IMAX = NI

              if (IKORD /= 0) then
                ! Restrict so (ij) <= (kl)
                IMAX = K
                JKINTOF = INTOF+(K-1)*NJ
                !XIJILS(1:NJ) = XIJKL(JKINTOF:JKINTOF+NJ-1)
                XIJILS(L:NL) = XIJKL(JKINTOF-1+L:JKINTOF-1+NL)
                XIJKL(JKINTOF-1+L) = Half*XIJKL(JKINTOF-1+L)
                XIJKL(JKINTOF+L:JKINTOF-1+NL) = Zero
              end if
              if ((IMAX >= 1) .and. (IMAX <= SKICKJ_TINY_N_MAX) .and. (NJ >= 1) .and. (NJ <= SKICKJ_TINY_K_MAX)) then
                call SKICKJ_MATML7_NN_TINY(SKII(ISOFF),CKJJ(ICOFF),XIJKL(INTOF),NKA,IMAX,NJ,FACS,FACTOR)
              else
                call MATML7(SKII(ISOFF),CKJJ(ICOFF),XIJKL(INTOF),NKA,IMAX,NKA,NJ,NJ,IMAX,FACS,FACTOR,0)
              end if
              if (IKORD /= 0) then
                XIJKL(JKINTOF-1+L:JKINTOF-1+NL) = XIJILS(L:NL)
                !XIJ(JKINTOF:JKINTOF+NJ-1) = XIJILS(1:NJ)
              end if

            end if
          end do

        end if
      end do
    end if
  end do
  ! (end over loop over Kb strings)
else if (IROUTE == 2) then
  ! S(I,Ka,Ib) = S(I,Ka,Ib) + sum(j) (ij!kl) C(j,Ka,Jb)
  do KB=1,NKB
    ! Number of nonvanishing connections from KB
    LL = 0
    KK = 0
    do L=1,NL
      if (KBJB(KB,L) /= 0) LL = LL+1
    end do
    do K=1,NK
      if (KBIB(KB,K) /= 0) KK = KK+1
    end do

    if ((KK /= 0) .and. (LL /= 0)) then
      do K=1,NK
        IB = KBIB(KB,K)
        if (IB /= 0) then
          SGNK = XKBIB(KB,K)
          do L=1,NL
            JB = KBJB(KB,L)
            if (JB /= 0) then
              SGNL = XKBJB(KB,L)
              FACTOR = SGNK*SGNL
              ! We have now a IB and Jb string, let's do it
              ISOFF = (IB-1)*NI*NKA+1
              ICOFF = (JB-1)*NJ*NKA+1
              INTOF = ((L-1)*NK+K-1)*NI*NJ+1

              if (IKORD /= 0) then
                ! Restrict so (ji) <= (kl)
                IKINTOF = INTOF+(K-1)*NI
                XIJILS(1:NI) = XIJKL(IKINTOF:IKINTOF+NI-1)
                XIJKL(IKINTOF-1+L) = Half*XIJKL(IKINTOF-1+L)
                XIJKL(IKINTOF+L:IKINTOF-1+NL) = Zero
              end if

              call MATML7(SKII(ISOFF),XIJKL(INTOF),CKJJ(ICOFF),NI,NKA,NI,NJ,NJ,NKA,FACS,FACTOR,0)

              if (IKORD /= 0) XIJKL(IKINTOF:IKINTOF+NI-1) = XIJILS(1:NI)

            end if
          end do
        end if
      end do
    end if
  end do
  ! (end over loop over Kb strings)

else if (IROUTE == 1) then
  write(u6,*) ' Sorry route 1 has retired, March 1997'
  !stop 'SKICKJ:Invalid route=1'
  call SYSABENDMSG('lucia_util/skickj','Internal error','')
  !do KB=1,NKB
  !  ! Number of nonvanishing a+lb !Kb>
  !  LL = 0
  !  do L=1,NL
  !    if (KBJB(KB,L) /= 0) LL = LL+1
  !  end do
  !
  !  IKEFF = 0
  !  do K=1,NK
  !    IB = KBIB(KB,K)
  !    if (IB == 0) cycle
  !    SGNK = XKBIB(KB,K)
  !
  !    if (IKORD == 0) then
  !       LI = NI
  !       IMIN = 1
  !    else
  !       LI = NI-K+1
  !       IMIN = K
  !    end if
  !
  !    do I=IMIN,NI
  !      IKEFF = IKEFF+1
  !      IOFF = (IKEFF-1)*NJ*LL
  !      ! Offset for S(1,IB,i)
  !      IBOFF(IKEFF) = (I-1)*NIB+IB
  !      LEFF = 0
  !      do L=1,NL
  !        JB = KBJB(KB,L)
  !        if (JB == 0) cycle
  !        LEFF = LEFF+1
  !        SGNL = XKBJB(KB,L)
  !        if ((IKORD == 1) .and. (I == K)) then
  !          FACTOR = Half*SGNK*SGNL
  !        else
  !          FACTOR = SGNK*SGNL
  !        end if
  !        JL0 = (LEFF-1)*NJ
  !        JLIK0 = (K-1)*NJ*NL*NI+(I-1)*NJ*NL+(L-1)*NJ
  !        do J=1,NJ
  !          JL = JL0+J
  !          ! Offsets for C(1,JB,j)
  !          JBOFF(JL) = (J-1)*NJB+JB
  !          ! integral * signs in SCR(jl,ik)
  !          ! Integrals are stored as (j l i k)
  !          SCR((IKEFF-1)*NJ*LL+JL) = FACTOR*XIJKL(JLIK)
  !          SCR(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
  !        end do
  !      end do
  !    end do
  !  end do
  !
  !  call GSAXPY_LUCIA(SKII,CKJJ,SCR,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)
  !end do
end if
! End of IROUTE branchning

contains

subroutine SKICKJ_MATML7_NN_TINY(C,A,B,M,N,K,FACTORC,FACTORAB)

  integer(kind=iwp), intent(in) :: M, N, K
  real(kind=wp), intent(inout) :: C(M,*)
  real(kind=wp), intent(in) :: A(M,*), B(K,*), FACTORC, FACTORAB
  integer(kind=iwp) :: I, J
  real(kind=wp) :: B1, B2, B3, B4

  select case (K)
    case (1)
      do J=1,N
        B1 = FACTORAB*B(1,J)
        if (FACTORC == Zero) then
          do I=1,M
            C(I,J) = B1*A(I,1)
          end do
        else
          do I=1,M
            C(I,J) = FACTORC*C(I,J)+B1*A(I,1)
          end do
        end if
      end do
    case (2)
      do J=1,N
        B1 = FACTORAB*B(1,J)
        B2 = FACTORAB*B(2,J)
        if (FACTORC == Zero) then
          do I=1,M
            C(I,J) = B1*A(I,1)+B2*A(I,2)
          end do
        else
          do I=1,M
            C(I,J) = FACTORC*C(I,J)+B1*A(I,1)+B2*A(I,2)
          end do
        end if
      end do
    case (3)
      do J=1,N
        B1 = FACTORAB*B(1,J)
        B2 = FACTORAB*B(2,J)
        B3 = FACTORAB*B(3,J)
        if (FACTORC == Zero) then
          do I=1,M
            C(I,J) = B1*A(I,1)+B2*A(I,2)+B3*A(I,3)
          end do
        else
          do I=1,M
            C(I,J) = FACTORC*C(I,J)+B1*A(I,1)+B2*A(I,2)+B3*A(I,3)
          end do
        end if
      end do
    case (4)
      do J=1,N
        B1 = FACTORAB*B(1,J)
        B2 = FACTORAB*B(2,J)
        B3 = FACTORAB*B(3,J)
        B4 = FACTORAB*B(4,J)
        if (FACTORC == Zero) then
          do I=1,M
            C(I,J) = B1*A(I,1)+B2*A(I,2)+B3*A(I,3)+B4*A(I,4)
          end do
        else
          do I=1,M
            C(I,J) = FACTORC*C(I,J)+B1*A(I,1)+B2*A(I,2)+B3*A(I,3)+B4*A(I,4)
          end do
        end if
      end do
  end select

end subroutine SKICKJ_MATML7_NN_TINY

end subroutine SKICKJ
