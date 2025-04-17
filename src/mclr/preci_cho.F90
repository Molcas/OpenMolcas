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
! Copyright (C) 2014, Mickael G. Delcey                                *
!***********************************************************************

subroutine Preci_cho(jS,nd,rOut,nbaj,fockii,fockai,focki,focka,fock,Sgn,A_J,nScr,iAdr)
!***********************************************************************
!                                                                      *
!     This routine remplaces preciaa, preciba and precibb              *
!     in case the new Cholesky algorithm is used,                      *
!     that is if only (ii|ab) and (ia|ib) integrals were computed      *
!                                                                      *
!     The code should be slightly more efficient as the integrals      *
!     are only read once and not for each distinct element             *
!                                                                      *
!     Written by M.G. Delcey, November 2014                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA
use input_mclr, only: LuChoInt, nAsh, nBas, nIsh, nOrb, nSym
use Constants, only: One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: jS, nd, nbaj, nScr
real(kind=wp), intent(inout) :: rout(*)
real(kind=wp), intent(in) :: fockii, fockai, focki(nbaj,nbaj), focka(nbaj,nbaj), fock(nbaj,nbaj), Sgn
real(kind=wp), intent(out) :: A_J(nScr)
integer(kind=iwp), intent(inout) :: iAdr
integer(kind=iwp) :: i, iAC, iBC, ijk, ilB, ip, ip1, iSym, iVB, jA, jAA, jB, jBB, jC, jCC, jD, jjA, jjB, jjC, jjD, jVert, kB, nlB, &
                     nO, nTri, nVirt, nVirt2
real(kind=wp) :: AABB, ACBB, BCBB, Factor, rDens, rDens1, rDens2, rFock

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)
nO = nAsh(jS)+nIsh(jS)
jVert = nOrb(jS)-nAsh(jS)-nIsh(jS)
nvirt = nOrb(jS)-nIsh(jS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Integral contribution

do iJK=1,2 ! 1=coulomb, 2=exchange
  factor = -One
  if (ijK == 2) factor = Three

  ! Read integrals

  do isym=1,nsym
    nvirt2 = nBas(isym)-nIsh(isym)
    call DDAFILE(LuChoInt(1),2,A_J,nvirt2**2,iAdr)

    ! iaa

    do jC=1,nAsh(isym)
      jjC = JC+nA(iSym)
      do jD=1,nAsh(isym)
        jjD = JD+nA(iSym)
        ip1 = (jD-1)*nvirt2+jC
        aabb = A_J(ip1) ! abab if exchange
        do jA=1,nAsh(jS)
          jjA = jA+nA(jS)
          do jB=1,jA
            jjB = jB+nA(jS)
            i = nTri-iTri(nd-jA+1,nd-jB+1)+1
            if (iJK == 1) then
              rDens1 = Two*Sgn*G2t(iTri(iTri(jjC,jjD),iTri(jjB,jjA)))
            else
              rDens1 = Four*Sgn*G2t((iTri(iTri(jjB,jjD),iTri(jjC,jjA))))
            end if
            rout(i) = rout(i)+rDens1*aabb
          end do
        end do
      end do
    end do

    if (iSym == jS) then
      do jA=1,nAsh(jS)
        jAA = jA+nA(jS)
        do jB=1,jA
          jBB = jB+nA(jS)
          i = nTri-iTri(nd-jA+1,nd-jB+1)+1
          do jC=1,nAsh(jS)
            jCC = jC+nA(jS)
            iBC = (jC-1)*nvirt+jB
            BCbb = A_J(iBC) ! BbCb if exchange
            iAC = (jC-1)*nvirt+jA
            ACbb = A_J(iAC) ! AbCb if exchange
            rDens1 = -Sgn*G1t(iTri(jAA,jCC))
            rDens2 = -Sgn*G1t(iTri(jBB,jCC))
            if (jAA == jCC) rDens1 = rdens1+Sgn
            if (jBB == jCC) rDens2 = rdens2+Sgn
            rout(i) = rout(i)+Two*rdens1*factor*BCbb+Two*rdens2*factor*ACbb

          end do
        end do
      end do

      ! iba

      do jA=1,nAsh(jS)
        ip = nTri-iTri(nd-ja+1,jVert)
        do jB=1,nAsh(jS)
          rDens = -Sgn*G1t(iTri(jA+nA(jS),jB+nA(jS)))
          if (jA == jB) rDens = rdens+Sgn*Two

          ivB = (jB-1)*nvirt+nAsh(jS)
          rOut(ip+1:ip+jVert) = rOut(ip+1:ip+jVert)+Two*factor*rDens*A_J(ivB+1:ivB+jVert)
        end do
      end do

      ! ibb

      i = nTri-nTri_Elem(jVert)+1
      do kB=nAsh(jS),nvirt-1
        nlB = nvirt-kb
        ilB = kB+1+nvirt*kb
        rout(i:i+nlB-1) = rout(i:i+nlB-1)+Sgn*Four*factor*A_J(ilB:ilB+nlB*nvirt-1:nvirt)
        i = i+nlB
      end do
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Fock matrix contribution

rFock = Sgn*(Fockii+FockAi)
!MGD tmp
!if (.false.) then
do jA=1,nAsh(jS)

  ! iaa

  jAA = jA+nA(jS)
  jjA = jA+nIsh(js)
  do jB=1,JA
    jBB = jB+nA(jS)
    jjB = jB+nIsh(js)
    i = nTri-iTri(nd-jA+1,nd-jB+1)+1
    rDens = G1t(iTri(jbb,jAA))
    rout(i) = rout(i)+Sgn*(Two*rdens*Fockii+Two*(Two*Focki(jjA,jjB)+Two*FockA(jjA,jjB)-Fock(jjB,jjA)))
  end do
  rout(i) = rout(i)-Four*rFock

  ! iba

  ip = nTri-iTri(nd-ja+1,nd-nAsh(js))
  rout(ip+1:ip+jVert) = rout(ip+1:ip+jVert)+Sgn*(Four*(FockI(nO+1:nO+jVert,ja+nIsh(js))+FockA(nO+1:nO+jVert,ja+nIsh(js)))- &
                                                 Fock(nO+1:nO+jVert,ja+nIsh(js)))
end do

! ibb

i = nTri-nTri_Elem(jVert)
do kB=nIsh(jS)+nAsh(jS),nOrb(jS)-1
  rOut(i+1) = rout(i+1)-Four*rFock
  rOut(i+1:i+nOrb(jS)-kB) = rout(i+1:i+nOrb(jS)-kB)+Sgn*Four*(Focki(kB+1,kB+1:nOrb(jS))+Focka(kB+1,kB+1:nOrb(jS)))
  i = i+nOrb(jS)-kB
end do
!end if

end subroutine Preci_cho
