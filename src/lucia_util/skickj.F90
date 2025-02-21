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
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SKICKJ(SKII,CKJJ,NKA,NKB,XIJKL,NI,NJ,NK,NL,MAXK,KBIB,XKBIB,KBJB,XKBJB,IKORD,FACS,IROUTE)
! Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
!          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
!
! Jeppe Olsen, Spring of 94
!
! : Note : Route 1 has retired, March 97

use lucia_data, only: MXPTSOB
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: SKII(*), XIJKL(*)
integer(kind=iwp), intent(in) :: NKA, NKB, NI, NJ, NK, NL, MAXK, KBIB(MAXK,*), KBJB(MAXK,*), IKORD, IROUTE
real(kind=wp), intent(in) :: CKJJ(*), XKBIB(MAXK,*), XKBJB(MAXK,*), FACS
integer(kind=iwp) :: IB, ICOFF, IKINTOF, IMAX, INTOF, ISOFF, JB, JKINTOF, K, KB, KK, L, LL
real(kind=wp) :: FACTOR, SGNK, SGNL, XIJILS(MXPTSOB)

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
              call MATML7(SKII(ISOFF),CKJJ(ICOFF),XIJKL(INTOF),NKA,IMAX,NKA,NJ,NJ,IMAX,FACS,FACTOR,0)
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

end subroutine SKICKJ
