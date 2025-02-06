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

implicit none
! Input
integer(kind=iwp) :: NKA, NKB, NI, NJ, NK, NL, MAXK, KBIB(MAXK,*), KBJB(MAXK,*), IKORD
real(kind=wp) :: SKII(*), CKJJ(*), RHO2B(*), XKBIB(MAXK,*), XKBJB(MAXK,*)
integer(kind=iwp) :: IB, ICOFF, IMAX, ISOFF, JB, K, KB, KK, KLOFF, L, LL
real(kind=wp) :: FACTOR, SGNK, SGNL

if (IKORD /= 0) then
  write(u6,*) ' ABTOR2 : IKORD /= 0'
  write(u6,*) ' I am not ready for this'
  !stop ' ABTOR2 : IKORD /= 0'
  call SYSABENDMSG('lucia_util/abtor2_gas','Internal error','')
end if

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
            !  do J=L,NL
            !    XIJILS(J) = XIJKL(JKINTOF-1+J)
            !  end do
            !  XIJKL(JKINTOF-1+L) = Half*XIJKL(JKINTOF-1+L)
            !  do J=L+1,NL
            !    XIJKL(JKINTOF-1+J) = Zero
            !  end do
            !end if
            call MATML7(RHO2B(KLOFF),SKII(ISOFF),CKJJ(ICOFF),NI,NJ,NKA,IMAX,NKA,NJ,One,FACTOR,1)
            !if (IKORD /= 0) then
            !  do J=L,NL
            !    XIJKL(JKINTOF-1+J) = XIJILS(J)
            !  end do
            !end if

          end if
        end do

      end if
    end do
  end if
end do
! (end over loop over Kb strings )

end subroutine ABTOR2
