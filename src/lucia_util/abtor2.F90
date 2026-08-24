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

subroutine ABTOR2(SKII,CKJJ,NKA,NKB,RHO2B,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD)
! Obtain contributions alpha-beta contributions to two-particle
! density matrix
!
! Rho2b(ij,kl)  = RHo2b(ij,kl)
!               + sum(Ka) Skii(Ka,i,Ib)<Ib!Eb(kl)!Jb> Ckjj(Ka,j,Jb)
!
! Jeppe Olsen, Fall of 96

use Constants, only: One
use Definitions, only: wp, iwp, u6
use lucia_parameters, only: ABTOR2_TINY_M_MAX, ABTOR2_TINY_N_MAX
#ifdef _CUDA_BLAS_
use, intrinsic :: iso_c_binding, only: c_int64_t
use ABTOR2_CUDA_INTERFACE, only: LUCIA_ABTOR2_CUDA_ROUTE
#endif

implicit none
integer(kind=iwp), intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, KBIB(MAXK,NK), KBJB(MAXK,NL), IKORD
real(kind=wp), intent(in) :: SKII(*), CKJJ(*), XKBIB(MAXK,NK), XKBJB(MAXK,NL)
real(kind=wp), intent(inout) :: RHO2B(NI*NJ*NK*NL)
integer(kind=iwp) :: IB, ICOFF, IMAX, ISOFF, JB, K, KB, KK, KLOFF, L, LL
real(kind=wp) :: FACTOR, SGNK, SGNL
#ifdef _CUDA_BLAS_
integer(c_int64_t) :: CudaStatus, NIBCUDA, NJBCUDA
#endif
if (IKORD /= 0) then
  write(u6,*) ' ABTOR2 : IKORD /= 0'
  write(u6,*) ' I am not ready for this'
  !stop ' ABTOR2 : IKORD /= 0'
  call SYSABENDMSG('lucia_util/abtor2_gas','Internal error','')
end if

#ifdef _CUDA_BLAS_
if ((NKA > 0) .and. (NKB > 0) .and. (NI > 0) .and. (NJ > 0) .and. (NK > 0) .and. &
    (NL > 0) .and. (MAXK > 0) .and. (NKB <= MAXK) .and. (IKORD == 0)) then
  NIBCUDA = maxval(KBIB(1:NKB,1:NK))
  NJBCUDA = maxval(KBJB(1:NKB,1:NL))
  if ((NIBCUDA > 0) .and. (NJBCUDA > 0)) then
    CudaStatus = LUCIA_ABTOR2_CUDA_ROUTE(RHO2B,SKII,CKJJ,NKA,NKB,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD, &
                                         NIBCUDA,NJBCUDA)
    if (CudaStatus == 1_c_int64_t) return
  end if
end if
#endif

! Excitations <Ib!Eb(kl)!Jb>
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
            KLOFF = ((L-1)*NK+K-1)*NI*NJ+1
            IMAX = NI

            !if (IKORD /= 0) then
            !  ! Restrict so (ij) <= (kl)
            !  IMAX = K
            !  JKINTOF = INTOF+(K-1)*NJ
            !  XIJILS(L:NL) = XIJKL(JKINTOF:JKINTOF+NL-1)
            !  XIJKL(JKINTOF-1+L) = Half*XIJKL(JKINTOF-1+L)
            !  XIJKL(JKINTOF+L:JKINTOF+NL-1) = Zero
            !end if
            if ((IMAX == NI) .and. (NI >= 1) .and. (NI <= ABTOR2_TINY_M_MAX) .and. &
                (NJ >= 1) .and. (NJ <= ABTOR2_TINY_N_MAX)) then
              call ABTOR2_MATML7_TN_TINY(RHO2B(KLOFF),SKII(ISOFF),CKJJ(ICOFF),NI,NJ,NKA,FACTOR)
            else
              call MATML7(RHO2B(KLOFF),SKII(ISOFF),CKJJ(ICOFF),NI,NJ,NKA,IMAX,NKA,NJ,One,FACTOR,1)
            end if
            !if (IKORD /= 0) XIJKL(JKINTOF+L-1:JKINTOF+NL-1) = XIJILS(L:NL)

          end if
        end do

      end if
    end do
  end if
end do
! (end over loop over Kb strings)

contains

subroutine ABTOR2_MATML7_TN_TINY(C,A,B,M,N,K,FACTORAB)
integer(kind=iwp), intent(in) :: M, N, K
real(kind=wp), intent(inout) :: C(M,*)
real(kind=wp), intent(in) :: A(K,*), B(K,*), FACTORAB
integer(kind=iwp) :: KA
real(kind=wp) :: B1, B2, B3, B4
real(kind=wp) :: S11, S21, S31, S41, S12, S22, S32, S42
real(kind=wp) :: S13, S23, S33, S43, S14, S24, S34, S44

S11 = 0.0_wp
S21 = 0.0_wp
S31 = 0.0_wp
S41 = 0.0_wp
S12 = 0.0_wp
S22 = 0.0_wp
S32 = 0.0_wp
S42 = 0.0_wp
S13 = 0.0_wp
S23 = 0.0_wp
S33 = 0.0_wp
S43 = 0.0_wp
S14 = 0.0_wp
S24 = 0.0_wp
S34 = 0.0_wp
S44 = 0.0_wp

select case (N)
case (1)
  do KA=1,K
    B1 = B(KA,1)
    S11 = S11+A(KA,1)*B1
    if (M >= 2) S21 = S21+A(KA,2)*B1
    if (M >= 3) S31 = S31+A(KA,3)*B1
    if (M >= 4) S41 = S41+A(KA,4)*B1
  end do
case (2)
  do KA=1,K
    B1 = B(KA,1)
    B2 = B(KA,2)
    S11 = S11+A(KA,1)*B1
    S12 = S12+A(KA,1)*B2
    if (M >= 2) then
      S21 = S21+A(KA,2)*B1
      S22 = S22+A(KA,2)*B2
    end if
    if (M >= 3) then
      S31 = S31+A(KA,3)*B1
      S32 = S32+A(KA,3)*B2
    end if
    if (M >= 4) then
      S41 = S41+A(KA,4)*B1
      S42 = S42+A(KA,4)*B2
    end if
  end do
case (3)
  do KA=1,K
    B1 = B(KA,1)
    B2 = B(KA,2)
    B3 = B(KA,3)
    S11 = S11+A(KA,1)*B1
    S12 = S12+A(KA,1)*B2
    S13 = S13+A(KA,1)*B3
    if (M >= 2) then
      S21 = S21+A(KA,2)*B1
      S22 = S22+A(KA,2)*B2
      S23 = S23+A(KA,2)*B3
    end if
    if (M >= 3) then
      S31 = S31+A(KA,3)*B1
      S32 = S32+A(KA,3)*B2
      S33 = S33+A(KA,3)*B3
    end if
    if (M >= 4) then
      S41 = S41+A(KA,4)*B1
      S42 = S42+A(KA,4)*B2
      S43 = S43+A(KA,4)*B3
    end if
  end do
case (4)
  do KA=1,K
    B1 = B(KA,1)
    B2 = B(KA,2)
    B3 = B(KA,3)
    B4 = B(KA,4)
    S11 = S11+A(KA,1)*B1
    S12 = S12+A(KA,1)*B2
    S13 = S13+A(KA,1)*B3
    S14 = S14+A(KA,1)*B4
    if (M >= 2) then
      S21 = S21+A(KA,2)*B1
      S22 = S22+A(KA,2)*B2
      S23 = S23+A(KA,2)*B3
      S24 = S24+A(KA,2)*B4
    end if
    if (M >= 3) then
      S31 = S31+A(KA,3)*B1
      S32 = S32+A(KA,3)*B2
      S33 = S33+A(KA,3)*B3
      S34 = S34+A(KA,3)*B4
    end if
    if (M >= 4) then
      S41 = S41+A(KA,4)*B1
      S42 = S42+A(KA,4)*B2
      S43 = S43+A(KA,4)*B3
      S44 = S44+A(KA,4)*B4
    end if
  end do
end select

C(1,1) = C(1,1)+FACTORAB*S11
if (M >= 2) C(2,1) = C(2,1)+FACTORAB*S21
if (M >= 3) C(3,1) = C(3,1)+FACTORAB*S31
if (M >= 4) C(4,1) = C(4,1)+FACTORAB*S41
if (N >= 2) then
  C(1,2) = C(1,2)+FACTORAB*S12
  if (M >= 2) C(2,2) = C(2,2)+FACTORAB*S22
  if (M >= 3) C(3,2) = C(3,2)+FACTORAB*S32
  if (M >= 4) C(4,2) = C(4,2)+FACTORAB*S42
end if
if (N >= 3) then
  C(1,3) = C(1,3)+FACTORAB*S13
  if (M >= 2) C(2,3) = C(2,3)+FACTORAB*S23
  if (M >= 3) C(3,3) = C(3,3)+FACTORAB*S33
  if (M >= 4) C(4,3) = C(4,3)+FACTORAB*S43
end if
if (N >= 4) then
  C(1,4) = C(1,4)+FACTORAB*S14
  if (M >= 2) C(2,4) = C(2,4)+FACTORAB*S24
  if (M >= 3) C(3,4) = C(3,4)+FACTORAB*S34
  if (M >= 4) C(4,4) = C(4,4)+FACTORAB*S44
end if
end subroutine ABTOR2_MATML7_TN_TINY

end subroutine ABTOR2
