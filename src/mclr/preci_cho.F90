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

subroutine Preci_cho(jS,nd,rOut,nbaj,fockii,fockai,focki,focka,fock,sign,A_J,nScr,iAdr)
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

use MCLR_Data, only: G1t, G2t
use MCLR_Data, only: nA
use input_mclr, only: nSym, nAsh, nIsh, nBas, nOrb, LuChoInt
use Constants, only: One, Two, Three, Four

implicit none
integer jS, nd
real*8 rout(*)
integer nbaj
real*8 fockii, fockai
real*8 focki(nbaj,nbaj), fock(nbaj,nbaj), focka(nbaj,nbaj)
real*8 sign
integer nScr
real*8 A_J(nScr)
integer iAdr
integer nTri, nO, jVert, nVirt, i1, ijk, iSym, nVirt2, jC, jjC, jD, jjD, ip1, jA, jjA, jB, jBB, jAA, iBC, iAC, ip, iVB, kB, nlB, &
        ilB, lB, jjB, jCC
real*8 Factor, AABB, rDens1, BCBB, ACBB, rDens2, rDens, rFock
! Statement functions
integer i, j, itri, itri1
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
itri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)
nO = nAsh(jS)+nIsh(jS)
jVert = nOrb(jS)-nAsh(jS)-nIsh(jS)
nvirt = nOrb(jS)-nIsh(jS)
i1 = nD-jVert+1
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
            i = itri1(jA,jB)
            if (iJK == 1) then
              rDens1 = Two*sign*G2t(itri(itri(jjC,jjD),itri(jjB,jjA)))
            else
              rDens1 = Four*sign*G2t((itri(itri(jjB,jjD),itri(jjC,jjA))))
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
          i = itri1(jA,jB)
          do jC=1,nAsh(jS)
            jCC = jC+nA(jS)
            iBC = (jC-1)*nvirt+jB
            BCbb = A_J(iBC) ! BbCb if exchange
            iAC = (jC-1)*nvirt+jA
            ACbb = A_J(iAC) ! AbCb if exchange
            rDens1 = -sign*G1t(itri(jAA,jCC))
            rDens2 = -sign*G1t(itri(jBB,jCC))
            if (jAA == jCC) rDens1 = rdens1+sign
            if (jBB == jCC) rDens2 = rdens2+sign
            rout(i) = rout(i)+Two*rdens1*factor*BCbb+Two*rdens2*factor*ACbb

          end do
        end do
      end do

      ! iba

      do jA=1,nAsh(jS)
        ip = itri1(ja,nd-jVert+1)
        do jB=1,nAsh(jS)
          rDens = -sign*G1t(iTri(jA+nA(jS),jB+nA(jS)))
          if (jA == jB) rDens = rdens+sign*Two

          ivB = (jB-1)*nvirt+nAsh(jS)+1
          call DaXpY_(jVert,Two*factor*rDens,A_J(ivB),1,rOut(ip),1)
        end do
      end do

      ! ibb

      i = itri1(i1,i1)
      do kB=nAsh(jS),nvirt-1
        nlB = nvirt-kb
        ilB = kB+1+nvirt*kb
        call daxpy_(nlB,sign*Four*factor,A_J(ilB),nvirt,rout(i),1)
        i = i+nlB
      end do
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Fock matrix contribution

rFock = sign*(Fockii+FockAi)
!MGD tmp
!if (.false.) then
do jA=1,nAsh(jS)

  ! iaa

  jAA = jA+nA(jS)
  jjA = jA+nIsh(js)
  do jB=1,JA
    jBB = jB+nA(jS)
    jjB = jB+nIsh(js)
    i = itri1(jA,jB)
    rDens = G1t(itri(jbb,jAA))
    rout(i) = rout(i)+Sign*(Two*rdens*Fockii+Two*(Two*Focki(jjA,jjB)+Two*FockA(jjA,jjB)-Fock(jjB,jjA)))
  end do
  rout(i) = rout(i)-Four*rFock

  ! iba

  ip = iTri1(ja,nAsh(js)+1)
  call DaXpY_(jVert,sign*Four,Focki(nO+1,ja+nIsh(js)),1,rout(ip),1)
  call DaXpY_(jVert,sign*Four,FockA(nO+1,ja+nIsh(js)),1,rout(ip),1)
  call DaXpY_(jVert,-sign,Fock(nO+1,ja+nIsh(js)),1,rout(ip),1)
end do

! ibb

i = itri1(i1,i1)-1
do kB=nIsh(jS)+nAsh(jS),nOrb(jS)-1
  rOut(i+1) = rout(i+1)-Four*rFock
  do lB=kb,nOrb(JS)-1
    i = i+1
    rOut(i) = rout(i)+sign*Four*Focki(kb+1,lb+1)+sign*Four*Focka(kb+1,lb+1)
  end do
end do
!end if

end subroutine Preci_cho
