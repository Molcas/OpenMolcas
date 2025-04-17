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

subroutine Preca_cho(iB,is,js,nd,rOut,nbaj,fockii,fockai,fockti,focki,focka,Sgn,A_J,nScr)
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

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA
use input_mclr, only: LuChoInt, nAsh, nBas, nIsh, nOrb, nSym
use Constants, only: One, Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iB, is, js, nd, nbaj, nScr
real(kind=wp), intent(inout) :: rout(nTri_Elem(nd))
real(kind=wp), intent(in) :: fockii, fockai, fockti, Focki(nbaj,nbaj), FockA(nBaj,nBaj), Sgn
real(kind=wp), intent(out) :: A_J(nScr)
integer(kind=iwp) :: i, iAdr, iBB, ii, ijk, ip, itAA, iU, iV, jA, jB, jCC, jDD, jVert, kSym, lSym, nI, nL, nO, nTri
real(kind=wp) :: Fact, Factor, Factor2, rDens, rDens1, rDens2, rDensabb, rDensabi, rDensabil, rDensabiu, rDensaii, rDensaiil, &
                 rDensaiiu, rf, rFock

!                                                                      *
!***********************************************************************
!                                                                      *
nO = nAsh(js)+nIsh(js)
nTri = nTri_Elem(nd)
iBB = ib+nA(is)
jVert = nOrb(js)-nIsh(js)-nAsh(js)
itAA = nTri_Elem(iBB)
!                                                                      *
!***********************************************************************
!                                                                      *
! Integral contribution
!
! iAdr=iAdr2 ! do not update iAdr2!
iAdr = 0
do ksym=1,nsym
  !ksym = Mul(jsym,)
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
              rDens2 = Two*Sgn*G2t(iTri(iTri(jCC,jDD),nTri_Elem(iBB)))
            else
              rDens2 = Two*Sgn*G2t(iTri(iTri(iBB,jDD),iTri(jCC,iBB)))
            end if

            rDensaii = factor*rDens2
            rDensabi = -factor2*rDens2
            rDensabb = factor2*rDens2
            rDensaiil = 0
            rDensaiiu = 0
            rDensabil = 0
            rDensabiu = 0

            if (iBB == jDD) then
              rDens1 = -Two*G1t(iTri(jCC,iBB))
              if (jCC == iBB) rdensaii = rdensaii-Two*factor
              rDensaiil = rDensaiil-factor*rDens1
              rDensabil = rDensabil+Sgn*rDens1
            else if ((iBB == jCC) .and. (iJK == 2)) then
              rDens1 = -Two*G1t(iTri(jDD,iBB))
              rDensaiiu = rDensaiiu-factor*rDens1
              rDensabiu = rDensabiu+Sgn*rDens1
            end if
            if ((iBB == jCC) .and. (iJK == 2)) then
              rDensaiil = rDensaiil-factor*14.0_wp*Sgn*G1t(iTri(jDD,iBB))
              rDensabil = rDensabil+Eight*Sgn*G1t(iTri(jDD,iBB))
              if (jDD == iBB) rdensaii = rdensaii+factor*14.0_wp*Sgn
            else if ((iBB == jDD) .and. (iJK == 2)) then
              rDensaiiu = rDensaiiu-factor*14.0_wp*Sgn*G1t(iTri(jCC,iBB))
              rDensabiu = rDensabiu+Eight*Sgn*G1t(iTri(jCC,iBB))
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
              ip = nTri-nTri_Elem(nd-ii+1)+1
              rout(ip:ip+ni-1) = rout(ip:ip+ni-1)+rDensaii*A_J((ii-1)*nl+ii:(ii-1)*nl+ii+ni-1)
              if ((iJK == 2) .and. (jCC > jDD)) rout(ip:ip+ni-1) = rout(ip:ip+ni-1)+rDensaiiu*A_J((ii-1)*nl+ii:(ii+ni-1)*nl+ii-1:nl)

              ! abi

              ip = nTri-iTri(nd-ii+1,jVert)+1
              rout(ip:ip+jvert-1) = rout(ip:ip+jvert-1)+rDensabi*A_J((ii-1)*nl+no+1:(ii-1)*nl+no+jvert)
              if ((iJK == 2) .and. (jCC > jDD)) &
                rout(ip:ip+jvert-1) = rout(ip:ip+jvert-1)+rDensabiu*A_J(no*nl+ii:(no+jvert)*nl+ii-1:nl)
            end do

            ! abb

            do ii=1,jvert
              ni = jvert-ii+1
              ip = nTri-nTri_Elem(nd-(nIsh(jS)+ii)+1)+1
              rout(ip:ip+ni-1) = rout(ip:ip+ni-1)+rDensabb*A_J((nO+ii-1)*nl+no+ii:(nO+ii-1)*nl+no+jvert)
              if ((iJK == 2) .and. (jCC > jDD)) &
                rout(ip:ip+ni-1) = rout(ip:ip+ni-1)+rDensabb*A_J((no+ii-1)*nl+no+ii:(no+ii+ni-1)*nl+no+ii-1:nl)
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

rFock = Sgn*(Two*(Fockii+Fockai)-Fockti)
rdens = Sgn*Two*G1t(nTri_Elem(ibb))

! aii

do jA=1,nIsh(jS)
  do jB=1,jA
    i = nTri-iTri(nd-ja+1,nd-jb+1)+1
    rout(i) = rout(i)-Sgn*Four*(Focka(jA,jB)+Focki(jA,jB))+rdens*Focki(ja,jb)
  end do
  rout(i) = rout(i)+Two*rfock

  ! abi

  ip = nTri-iTri(nd-jA+1,jVert)
  Fact = (Two-Two*G1t(itAA))
  rOut(ip+1:ip+jVert) = rOut(ip+1:ip+jVert)+Sgn*(Fact*FockI(nO+1:nO+jVert,jA)+Two*FockA(nO+1:nO+jVert,jA))
end do

! abb

ip = nTri-nTri_Elem(jVert)+1
rF = Sgn*Fockti
do iI=nAsh(js)+nIsh(js)+1,nBas(js)
  rOut(ip) = rout(ip)-Two*rF+rDens*FockI(iI,ii)
  ip = ip+1
  rOut(ip:ip+nBas(jS)-iI-1) = rOut(ip:ip+nBas(jS)-iI-1)+rDens*FockI(iI,iI+1:nBas(js))
  ip = ip+nBas(js)-iI
end do

end subroutine Preca_cho
