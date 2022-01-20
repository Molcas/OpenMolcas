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

function VExch(ZP,NP,ZQ,NQ,LA,nProj,iCoShll)
!***********************************************************************
!                                                                      *
!     VExch calculates the atomic integral                             *
!     <zp,np,la | Sum(core) Kc |zq,nq,la>                              *
!                                                                      *
!***********************************************************************

use Basis_Info, only: Shells
use AMatrix, only: DFAC, KOSUU, NYU, RCA
use Constants, only: Zero, One, Two, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: VExch
real(kind=wp), intent(in) :: ZP, ZQ
integer(kind=iwp), intent(in) :: NP, NQ, LA, nProj, iCoShll
integer(kind=iwp) :: ICORB, iCoSh, INU, ISIM, IT1, IT2, IT3, IT4, K, KOMAX, L, L1, L2, LMT, nExpon, NR, NS, NU, NUT
real(kind=wp) :: CR, CS, DL2, FOcc, OrbPS, RCAT, RTT1, RTT2, RTT3, RTT4, RTT5, RTT6, RTT7, SUMA, VPQ, VR, VS, ZR, ZS
real(kind=wp), parameter :: SQRT2PI = sqrt(Two/Pi)

if (nProj > 4) then
  write(u6,*) 'VExch: nProj',nProj
  write(u6,*) 'Abend: Implementation ready only up to g-core.'
  write(u6,*) '       Update common block /CONST/.'
  call Abend()
end if
if ((NP > 5) .or. (NQ > 5)) then
  write(u6,*) 'VExch: NP,NQ',NP,NQ
  write(u6,*) 'Abend: Implementation ready only up to g-valence.'
  write(u6,*) '       Update common block /CONST/.'
  call Abend()
end if

!PIPPI = (Half/PI)**0.75_wp
L1 = LA
VPQ = VF(2*NP,ZP)*VF(2*NQ,ZQ)
!PSMT = Zero
!ls start adding
Vexch = Zero
!ls end adding
! loop over angular momentum
do ISIM=1,nProj+1
  iCoSh = iCoShll+ISIM-1
  !EDUM = Half*real(ISIM,kind=wp)+Quart
  NR = ISIM
  NS = ISIM
  L2 = ISIM
  DL2 = real(2*(L2-1)+1,kind=wp)
  LMT = ((L1-1)*L1)/2+L2
  if (L1 < L2) LMT = ((L2-1)*L2)/2+L1
  KOMAX = KOSUU(LMT)
  ! loop over core orbitals of a given angular momentum
  do ICORB=1,Shells(iCoSh)%nBasis
    OrbPS = Zero
    do INU=1,KOMAX
      NUT = NYU(INU,LMT)
      NU = NUT
      RCAT = RCA(INU,LMT)*DL2
      !      BEGIN KSM.
      ! COMPUTATION OF K SUPER MATRIX OVER CGTO
      ! NU IS UPPER INDEX NYU IN ROOTHAAN'S NOTATION

      IT1 = NP+NR-NU-1
      IT2 = NQ+NS+NU
      IT3 = NQ+NS-NU-1
      IT4 = NP+NR+NU
      SUMA = Zero
      nExpon = Shells(iCoSh)%nExp
      do K=1,nExpon
        ZR = Shells(iCoSh)%Exp(K)
        CR = Shells(iCoSh)%Cff_c(K,ICORB,2)
        !LS coeff of unnormalized non-diagonal cartesian GTF
        ! to be used in ecpaimp       CR = CR/(PIPPI*(Four*ZR)**EDUM)
        VR = VF(2*NR,ZR)
        RTT1 = Half*(ZP+ZR)
        do L=1,nExpon
          ZS = Shells(iCoSh)%Exp(L)
          CS = Shells(iCoSh)%Cff_c(L,ICORB,2)
          ! to be used in ecpaimp          CS = CS/(PIPPI*(Four*ZS)**EDUM)
          VS = VF(2*NS,ZS)
          RTT2 = Half*(ZQ+ZS)
          RTT3 = RTT1/RTT2
          RTT4 = One/RTT3
          call AUXC((IT1+1)/2,IT2,RTT3,RTT5)
          call AUXC((IT3+1)/2,IT4,RTT4,RTT6)
          RTT7 = VF(IT1,RTT1)*VF(IT2,RTT2)*RTT5+VF(IT3,RTT2)*VF(IT4,RTT1)*RTT6
          SUMA = SUMA+RTT7*CR*CS/sqrt(VPQ*VR*VS)
        end do
      end do
      !      END KSM.
      !ls PSMT = PSMT+RCAT*SQRT2PI*SUMA
      !ls start adding
      OrbPS = OrbPS+RCAT*SQRT2PI*SUMA
      !ls end adding
    end do
    !ls start adding
    FOcc = Shells(iCoSh)%Occ(iCOrb)
    Vexch = Vexch+Two*OrbPS*FOcc
    !ls start adding
  end do ! end of loop over core orbitals of a given angular momentum
end do ! end of loop over angular momentum
!ls VEXCH = Two*PSMT

return

contains

function VF(N,X)

  real(kind=wp) :: VF
  integer(kind=iwp), intent(in) :: N
  real(kind=wp), intent(in) :: X

  VF = DFAC(N)/sqrt(X)**(N+1)

end function VF

end function VExch
