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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine SKICKJ_MCLR(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,IROUTE)
! Calculate S(Ka,Ib,i) = S(Ka,Ib,i)
!          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
!
! Jeppe Olsen, Spring of 94

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NKA, NIB, NJB, NKB, NI, NJ, NK, NL, MAXK, KBIB(MAXK,*), KBJB(MAXK,*), IKORD, IROUTE
real(kind=wp), intent(inout) :: SKII(NKA*NI,*)
real(kind=wp), intent(in) :: CKJJ(NKA*NJ,*), XIJKL(*), XKBIB(MAXK,*), XKBJB(MAXK,*)
integer(kind=iwp) :: I, IB, IKEFF, IMIN, INTOF, IOFF, J, JB, JL, JL0, JLIK0, K, KB, KK, L, LEFF, LENGTH, LL, MAXORB
real(kind=wp) :: FACTOR, SGNK, SGNL
integer(kind=iwp), allocatable :: IBOFF(:), JBOFF(:)
real(kind=wp), allocatable :: KSKICK(:)
integer(kind=iwp), parameter :: MXTSOB = 35

! Note if Iroute = 2 the form is S(i,Ka,Ib)
! Note if Iroute = 2 the form is C(j,Ka,Jb)

MAXORB = max(NI,NJ,NK,NL)
LENGTH = MAXORB*MAXORB*MAXORB*MAXORB
call mma_allocate(KSKICK,LENGTH,Label='KSKICK')
call mma_allocate(IBOFF,MXTSOB*MXTSOB,Label='IBOFF')
call mma_allocate(JBOFF,MXTSOB*MXTSOB,Label='JBOFF')

if (MAXORB > MXTSOB) then
  write(u6,*) ' SKICKJ_MCLR : Too many orbs : NI > MXTSOB'
  write(u6,*) ' NI, MXTSOB ',max(NI,NJ,NK,NL),MXTSOB
  write(u6,*) ' Redim MXTSOB in SKICKJ_MCLR'
  call Abend()
end if

if (IROUTE == 3) then
  ! S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
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
              !ISOFF = (IB-1)*NI*NKA+1
              !ICOFF = (JB-1)*NJ*NKA+1
              INTOF = ((L-1)*NK+K-1)*NI*NJ+1

              call DGEMM_('N','N',NKA,NI,NJ,FACTOR,CKJJ(1,jB),max(1,NKA),XIJKL(INTOF),max(1,NJ),ONE,SKII(1,iB),max(1,NKA))

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
              !ISOFF = (IB-1)*NI*NKA+1
              !ICOFF = (JB-1)*NJ*NKA+1
              INTOF = ((L-1)*NK+K-1)*NI*NJ+1

              call DGEMM_('N','N',NI,NKA,NJ,FACTOR,XIJKL(INTOF),max(1,NI),CKJJ(1,JB),max(1,NJ),ONE,SKII(1,IB),max(1,NI))
            end if
          end do
        end if
      end do
    end if
  end do
  ! (end over loop over Kb strings)

else if (IROUTE == 1) then

  do KB=1,NKB
    ! Number of nonvanishing a+lb !Kb>
    LL = 0
    do L=1,NL
      if (KBJB(KB,L) /= 0) LL = LL+1
    end do

    IKEFF = 0
    do K=1,NK
      IB = KBIB(KB,K)
      if (IB == 0) cycle
      SGNK = XKBIB(KB,K)

      if (IKORD == 0) then
        IMIN = 1
      else
        IMIN = K
      end if

      do I=IMIN,NI
        IKEFF = IKEFF+1
        IOFF = (IKEFF-1)*NJ*LL
        ! Offset for S(1,IB,i)
        IBOFF(IKEFF) = (I-1)*NIB+IB
        LEFF = 0
        do L=1,NL
          JB = KBJB(KB,L)
          if (JB == 0) cycle
          LEFF = LEFF+1
          SGNL = XKBJB(KB,L)
          if ((IKORD == 1) .and. (I == K)) then
            FACTOR = Half*SGNK*SGNL
          else
            FACTOR = SGNK*SGNL
          end if
          JL0 = (LEFF-1)*NJ
          JLIK0 = (K-1)*NJ*NL*NI+(I-1)*NJ*NL+(L-1)*NJ

          do J=1,NJ
            JL = JL0+J
            ! Offsets for C(1,JB,j)
            JBOFF(JL) = (J-1)*NJB+JB
            ! integral * signs in SCR(jl,ik)
            ! Integrals are stored as (j l i k)
            KSKICK(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
          end do

        end do

      end do

    end do

    call GSAXPY(SKII,CKJJ,KSKICK,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)

  end do
end if
! End of IROUTE branching

call mma_deallocate(KSKICK)
call mma_deallocate(IBOFF)
call mma_deallocate(JBOFF)

end subroutine SKICKJ_MCLR
