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

subroutine Preca_cho(iB,is,js,nd,ir,rOut,nbai,nbaj,fockii,fockai,fockti,focki,focka,fock,sign,A_J,A_K,Scr,nScr,iAdr)
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

implicit none
integer iB, is, js, nd, ir
real*8 rout(nd*(nd+1)/2)
integer nbai, nbaj
real*8 fockii, fockai, fockti
real*8 Fock(nbaj,nbaj), FockA(nBaj,nBaj), Focki(nbaj,nbaj)
real*8 Sign
integer nScr
real*8 A_J(nScr), A_K(nScr), Scr(nScr)
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
        factor = 2.0d0
        factor2 = 1.0d0
        if (iJK == 2) factor = 1.0d0 ! exchange
        if (iJK == 2) factor2 = 2.0d0 ! exchange
        do lsym=1,nsym
          nl = nBas(lsym)
          if (nl <= 0) cycle
          call DDAFILE(LuChoInt(2),2,A_J,nl**2,iAdr)

          if (lsym == js) then
            if (iJK == 1) then
              rDens2 = 2.0d0*sign*G2t(itri(itri(jCC,jDD),itri(iBB,iBB)))
            else
              rDens2 = 2.0d0*sign*G2t(itri(itri(iBB,jDD),itri(jCC,iBB)))
            end if

            rDensaii = factor*rDens2
            rDensabi = -factor2*rDens2
            rDensabb = factor2*rDens2
            rDensaiil = 0
            rDensaiiu = 0
            rDensabil = 0
            rDensabiu = 0

            if (iBB == jDD) then
              rDens1 = -2.0d0*G1t(itri(jCC,iBB))
              if (jCC == iBB) rdensaii = rdensaii-2.0d0*factor
              rDensaiil = rDensaiil-factor*rDens1
              rDensabil = rDensabil+sign*rDens1
            else if ((iBB == jCC) .and. (iJK == 2)) then
              rDens1 = -2.0d0*G1t(itri(jDD,iBB))
              rDensaiiu = rDensaiiu-factor*rDens1
              rDensabiu = rDensabiu+sign*rDens1
            end if
            if ((iBB == jCC) .and. (iJK == 2)) then
              rDensaiil = rDensaiil-factor*14.0d0*sign*G1t(itri(jDD,iBB))
              rDensabil = rDensabil+8.0d0*sign*G1t(itri(jDD,iBB))
              if (jDD == iBB) rdensaii = rdensaii+factor*14.0d0*sign
            else if ((iBB == jDD) .and. (iJK == 2)) then
              rDensaiiu = rDensaiiu-factor*14.0d0*sign*G1t(itri(jCC,iBB))
              rDensabiu = rDensabiu+8.0d0*sign*G1t(itri(jCC,iBB))
            end if
            rDensaiiu = rDensaiiu+rDensaii
            rDensaii = rDensaiil+rDensaii
            rDensabiu = rDensabiu+rDensabi
            rDensabi = rDensabil+rDensabi
            if ((iJK == 1) .and. (jCC > jDD)) then
              rDensaii = rDensaii*2.0d0
              rDensabi = rDensabi*2.0d0
              rDensabb = rDensabb*2.0d0
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

rFock = sign*2.0d0*Fockii+sign*2.0d0*Fockai-sign*Fockti
rdens = sign*2.0d0*G1t(itri(ibb,ibb))

! aii

do jA=1,nIsh(jS)
  do jB=1,jA
    i = itri1(ja,jb)
    rout(i) = rout(i)-sign*4.0d0*(Focka(jA,jB)+Focki(jA,jB))+rdens*Focki(ja,jb)
  end do
  rout(i) = rout(i)+2.0d0*rfock

  ! abi

  ip = itri1(jA,nd-jVert+1)
  Fact = (2.0d0-2.0d0*G1t(itAA))
  call DaxPy_(jVert,Sign*Fact,FockI(nO+1,jA),1,rOut(ip),1)
  Fact = 2.0d0
  call DaxPy_(jVert,Sign*Fact,FockA(nO+1,jA),1,rOut(ip),1)
end do

! abb

ip = iTri1(i2,i2)
rF = sign*Fockti
do iI=nAsh(js)+nIsh(js)+1,nBas(js)
  rOut(ip) = rout(ip)-2.0d0*rF+rDens*FockI(iI,ii)
  ip = ip+1
  do iJ=iI+1,Nbas(js)
    rOut(ip) = rout(ip)+rDens*FockI(iI,iJ)
    ip = ip+1
  end do
end do

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(ir)
  call Unused_integer(nbai)
  call Unused_real_array(fock)
  call Unused_real_array(A_K)
  call Unused_real_array(Scr)
end if

end subroutine Preca_cho
