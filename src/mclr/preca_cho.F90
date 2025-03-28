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

subroutine Preca_cho(iB,is,js,nd,rOut,nbaj,fockii,fockai,fockti,focki,focka,sign,A_J,nScr,iAdr)
!***********************************************************************
!                                                                      *
!     This routine replaces precaii, precabi and precabb               *
!     in case the new Cholesky alrgorithm is used,                     *
!     that is if only (ii|ab) and (ia|ib) integrals were computed      *
!                                                                      *
!     The code should be slightly more efficient as the integrals      *
!     are only read once and not for each distinct element             *
!                                                                      *
!     Written by M.G. Delcey, november 2014                            *
!                                                                      *
!***********************************************************************

use Arrays, only: G1t, G2t
use MCLR_Data, only: nA
use input_mclr, only: nSym, nAsh, nIsh, nBas, nOrb, LuChoInt
use Constants, only: One, Two, Four, Eight
use Definitions, only: wp

implicit none
integer iB, is, js, nd
real*8 rout(nd*(nd+1)/2)
integer nbaj
real*8 fockii, fockai, fockti
real*8 FockA(nBaj,nBaj), Focki(nbaj,nbaj)
real*8 Sign
integer nScr
real*8 A_J(nScr)
integer nTri, nO, iBB, jVert, itAA, i2, iAdr, kSym, iV, jCC, iU, jDD, ijk, lSym, nL, ii, nI, ip, jA, jB, ij
real*8 Factor, Factor2, rDens2, rDensaii, rDensabi, rDensabb, rf, rDens, rDensaiil, rDensaiiu, rDensabil, rDensabiu, rFock, &
       rDens1, Fact
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nO = nAsh(js)+nIsh(js)
nTri = itri(nd,nd)
iBB = ib+nA(is)
jVert = nOrb(js)-nIsh(js)-nAsh(js)
itAA = itri(iBB,iBB)
i2 = nD-jVert+1
!                                                                      *
!***********************************************************************
!                                                                      *
! Integral contribution
!
! iAdr=iAdr2 ! do not update iAdr2!
iAdr = 0
do ksym=1,nsym
  !ksym=iEor(jsym,)
  !!!!!!!!!!!!!
  do iv=1,nAsh(kSym)
    jCC = iv+nA(kSym)
    do iu=1,iv
      jDD = iu+nA(kSym)
      do iJK=1,2 ! 1=coulomb, 2=exchange
        factor = Two
        factor2 = One
        if (iJK == 2) factor = One  ! exchange
        if (iJK == 2) factor2 = Two ! exchange
        do lsym=1,nsym
          nl = nBas(lsym)
          if (nl <= 0) cycle
          call DDAFILE(LuChoInt(2),2,A_J,nl**2,iAdr)

          if (lsym == js) then
            if (iJK == 1) then
              rDens2 = Two*sign*G2t(itri(itri(jCC,jDD),itri(iBB,iBB)))
            else
              rDens2 = Two*sign*G2t(itri(itri(iBB,jDD),itri(jCC,iBB)))
            end if

            rDensaii = factor*rDens2
            rDensabi = -factor2*rDens2
            rDensabb = factor2*rDens2
            rDensaiil = 0
            rDensaiiu = 0
            rDensabil = 0
            rDensabiu = 0

            if (iBB == jDD) then
              rDens1 = -Two*G1t(itri(jCC,iBB))
              if (jCC == iBB) rdensaii = rdensaii-Two*factor
              rDensaiil = rDensaiil-factor*rDens1
              rDensabil = rDensabil+sign*rDens1
            else if ((iBB == jCC) .and. (iJK == 2)) then
              rDens1 = -Two*G1t(itri(jDD,iBB))
              rDensaiiu = rDensaiiu-factor*rDens1
              rDensabiu = rDensabiu+sign*rDens1
            end if
            if ((iBB == jCC) .and. (iJK == 2)) then
              rDensaiil = rDensaiil-factor*14.0_wp*sign*G1t(itri(jDD,iBB))
              rDensabil = rDensabil+Eight*sign*G1t(itri(jDD,iBB))
              if (jDD == iBB) rdensaii = rdensaii+factor*14.0_wp*sign
            else if ((iBB == jDD) .and. (iJK == 2)) then
              rDensaiiu = rDensaiiu-factor*14.0_wp*sign*G1t(itri(jCC,iBB))
              rDensabiu = rDensabiu+Eight*sign*G1t(itri(jCC,iBB))
            end if
            rDensaiiu = rDensaiiu+rDensaii
            rDensaii = rDensaiil+rDensaii
            rDensabiu = rDensabiu+rDensabi
            rDensabi = rDensabil+rDensabi
            if ((iJK == 1) .and. (jCC > jDD)) then
              rDensaii = rDensaii*Two
              rDensabi = rDensabi*Two
              rDensabb = rDensabb*Two
            end if

            do ii=1,nIsh(js)

              ! aii

              ni = nIsh(jS)-ii+1
              ip = iTri1(ii,ii)
              call DaXpY_(ni,rDensaii,A_J((ii-1)*nl+ii),1,rout(ip),1)
              if ((iJK == 2) .and. (jCC > jDD)) call DaXpY_(ni,rDensaiiu,A_J((ii-1)*nl+ii),nl,rout(ip),1)

              ! abi

              ip = itri1(ii,nd-jVert+1)
              call DaXpY_(jvert,rDensabi,A_J((ii-1)*nl+no+1),1,rout(ip),1)
              if ((iJK == 2) .and. (jCC > jDD)) call DaXpY_(jvert,rDensabiu,A_J(no*nl+ii),nl,rout(ip),1)
            end do

            ! abb

            do ii=1,jvert
              ni = jvert-ii+1
              ip = itri1(nIsh(jS)+ii,nIsh(jS)+ii)
              call DaXpY_(ni,rDensabb,A_J((nO+ii-1)*nl+no+ii),1,rout(ip),1)
              if ((iJK == 2) .and. (jCC > jDD)) call DaXpY_(ni,rDensabb,A_J((nO+ii-1)*nl+no+ii),nl,rout(ip),1)
            end do
          end if

        end do
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Fock matrix contribution

rFock = sign*Two*Fockii+sign*Two*Fockai-sign*Fockti
rdens = sign*Two*G1t(itri(ibb,ibb))

! aii

do jA=1,nIsh(jS)
  do jB=1,jA
    i = itri1(ja,jb)
    rout(i) = rout(i)-sign*Four*(Focka(jA,jB)+Focki(jA,jB))+rdens*Focki(ja,jb)
  end do
  rout(i) = rout(i)+Two*rfock

  ! abi

  ip = itri1(jA,nd-jVert+1)
  Fact = (Two-Two*G1t(itAA))
  call DaxPy_(jVert,Sign*Fact,FockI(nO+1,jA),1,rOut(ip),1)
  Fact = Two
  call DaxPy_(jVert,Sign*Fact,FockA(nO+1,jA),1,rOut(ip),1)
end do

! abb

ip = iTri1(i2,i2)
rF = sign*Fockti
do iI=nAsh(js)+nIsh(js)+1,nBas(js)
  rOut(ip) = rout(ip)-Two*rF+rDens*FockI(iI,ii)
  ip = ip+1
  do iJ=iI+1,Nbas(js)
    rOut(ip) = rout(ip)+rDens*FockI(iI,iJ)
    ip = ip+1
  end do
end do

end subroutine Preca_cho
