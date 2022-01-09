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

real*8 function VExch(ZP,NP,ZQ,NQ,LA,nProj,iCoShll)
!***********************************************************************
!                                                                      *
!     VExch calculates the atomic integral                             *
!     <zp,np,la | Sum(core) Kc |zq,nq,la>                              *
!                                                                      *
!***********************************************************************

use Basis_Info
implicit real*8(A-H,O-Z)
! auxiliary constant pool:       ready only up to g-valence/g-core
#include "const.fh"
! Define statement function
VF(N,X) = DFAC(N)/sqrt(X)**(N+1)

if (nProj > 4) then
  write(6,*) 'VExch: nProj',nProj
  write(6,*) 'Abend: Implementation ready only up to g-core.'
  write(6,*) '       Update common block /CONST/.'
  call Abend()
end if
if ((NP > 5) .or. (NQ > 5)) then
  write(6,*) 'VExch: NP,NQ',NP,NQ
  write(6,*) 'Abend: Implementation ready only up to g-valence.'
  write(6,*) '       Update common block /CONST/.'
  call Abend()
end if

!PI = 2.D0*ACOS(0.D0)
!PIPPI = (0.5D0/PI)**0.75D0
L1 = LA
VPQ = VF(2*NP,ZP)*VF(2*NQ,ZQ)
!PSMT = 0.D0
!ls start adding
Vexch = 0.d0
!ls end adding
! loop over angular momentum
do ISIM=1,nProj+1
  iCoSh = iCoShll+ISIM-1
  !EDUM = 0.5D0*DBLE(ISIM)+0.25D0
  NR = ISIM
  NS = ISIM
  L2 = ISIM
  DL2 = dble(2*(L2-1)+1)
  LMT = ((L1-1)*L1)/2+L2
  if (L1 < L2) LMT = ((L2-1)*L2)/2+L1
  KOMAX = KOSUU(LMT)
  ! loop over core orbitals of a given angular momentum
  do ICORB=1,Shells(iCoSh)%nBasis
    OrbPS = 0d0
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
      SUMA = 0.d0
      nExpon = Shells(iCoSh)%nExp
      do K=1,nExpon
        ZR = Shells(iCoSh)%exp(K)
        CR = Shells(iCoSh)%Cff_c(K,ICORB,2)
        !LS coeff of unnormalized non-diagonal cartesian GTF
        ! to be used in ecpaimp       CR=CR/(PIPPI*(4.D0*ZR)**EDUM)
        VR = VF(2*NR,ZR)
        RTT1 = 0.5d0*(ZP+ZR)
        do L=1,nExpon
          ZS = Shells(iCoSh)%exp(L)
          CS = Shells(iCoSh)%Cff_c(L,ICORB,2)
          ! to be used in ecpaimp          CS=CS/(PIPPI*(4.D0*ZS)**EDUM)
          VS = VF(2*NS,ZS)
          RTT2 = 0.5d0*(ZQ+ZS)
          RTT3 = RTT1/RTT2
          RTT4 = 1.d0/RTT3
          call AUXC((IT1+1)/2,IT2,RTT3,RTT5)
          call AUXC((IT3+1)/2,IT4,RTT4,RTT6)
          RTT7 = VF(IT1,RTT1)*VF(IT2,RTT2)*RTT5+VF(IT3,RTT2)*VF(IT4,RTT1)*RTT6
          SUMA = SUMA+RTT7*CR*CS/sqrt(VPQ*VR*VS)
        end do
      end do
      !      END KSM.
      !ls PSMT = PSMT+RCAT*0.797884561D0*SUMA
      !ls start adding
      OrbPS = OrbPS+RCAT*0.797884561d0*SUMA
      !ls end adding
    end do
    !ls start adding
    FOcc = Shells(iCoSh)%Occ(iCOrb)
    Vexch = Vexch+2d0*OrbPS*FOcc
    !ls start adding
  end do ! end of loop over core orbitals of a given angular momentum
end do ! end of loop over angular momentum
!ls VEXCH = 2.D0*PSMT

return

end function VExch
