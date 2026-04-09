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

subroutine SetUp_MCLR(DSYM)
! Setup pointers and size of matrices (includes in Pointers.fh)

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipCM, ipMat, ipMatBA, ipMatLT, ipMO, n1Dens, n2Dens, nA, nB, nCMO, nDens, nDensC, nMBA, nNA, pInt1, pInt2
use input_mclr, only: iMethod, nAsh, nBas, nDel, nFro, nIsh, nOrb, nRS1, nRS2, nRs3, nSym, PT2, TimeDep
use caspt2_module, only: nBMx
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: dsym
integer(kind=iwp) :: i, i1, iExt0, iExt1, iExt2, iExt3, iExt4, iInt0, iInt1, iInt2, iInt3, iInt4, iiSym, ijNum, ijOrb, ijS, ijSym, &
                     iOff, iOrb, ip, iPlus, ipP, iS, jjSym, jOrb, jS, kkSym, klNum, klOrb, klS, klSym, kOrb, kS, llSym, lOrb, lS, &
                     mATAB, nDensLT
integer(kind=iwp), external :: iPntSO

!                                                                      *
!***********************************************************************
!                                                                      *
ip = 1
ipMO(:,:,:) = 0

nBMx = max(0,maxval(nBas(1:nsym)))
nB(1:nsym) = nAsh(1:nsym)+nIsh(1:nsym)
nA(1) = 0
do iS=2,nsym
  nA(iS) = nA(iS-1)+nAsh(iS-1)
end do

nna = sum(nAsh(1:nsym))
n2Dens = nTri_Elem(nnA**2)

n1Dens = nnA*nnA
ip = 1
do kS=1,nSym
  do lS=1,nsym
    klS = Mul(lS,kS)
    do jS=1,nSym
      do iS=1,nSym
        ijS = Mul(iS,jS)
        ipp = nOrb(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
        if ((Mul(ijS,klS) == DSym) .and. (ipp /= 0)) then
          ipMO(kS,lS,jS) = ip
          ip = ip+ipp
        end if
      end do
    end do
  end do
end do

nmba = ip-1
nDens = 0
nDensLT = 0
nCMO = 0
nDensC = 0
ipMat(:,:) = 0
ipMatLT(:,:) = 0
ipCM(:) = 0

do jS=1,nSym
  do iS=1,js
    if (Mul(iS,jS) == DSym) then
      if (is < js) then
        ipMatLT(jS,iS) = ipntso(js-1,is-1,ibset(0,dsym-1),nBas)+1
        nDensLT = nDensLT+nBas(iS)*nBas(jS)
        iExt0 = nOrb(is)-nIsh(is)
        iExt1 = nOrb(iS)-nRs1(iS)
        iExt2 = nOrb(iS)-nRs2(iS)
        iExt3 = nOrb(iS)-nRs3(iS)
        iInt4 = nOrb(js)-nish(js)-nAsh(js)
        iExt4 = nIsh(is)+nAsh(is)
        nDensC = nDensC+iExt0*nIsh(js)+iExt1*nRs1(js)+iExt2*nRs2(js)+iExt3*nRs3(js)+iExt4*iInt4
      end if
      if (is == js) then
        i1 = nOrb(is)-nish(is)-nAsh(is)
        ipMatLT(jS,iS) = nDensLT+1
        nDensLT = nDensLT+nTri_Elem(nBas(iS))
        iint0 = nOrb(is)-nIsh(is)
        iint1 = i1+nRs2(is)+nRs3(is)
        iint2 = i1+nRs3(is)
        iint3 = i1
        nDensC = nDensC+iint0*nIsh(is)+iint1*nRs1(is)+iint2*nRs2(is)+iint3*nRs3(is)
      end if
    end if
  end do
  ipCM(jS) = nCMO+1
  nCMO = nCMO+nBas(js)**2
end do

if (TimeDep) nDensC = 2*nDensC
nDens = 0
matab = 1
do jS=1,nSym
  do iS=1,nSym
    if (Mul(iS,jS) == DSym) then
      ipMatba(is,js) = matab
      ipMat(jS,iS) = nDens+1
      nDens = nDens+nBas(iS)*nBas(jS)
      matab = matab+nash(js)*nOrb(iS)
    end if
  end do
end do

! To begin with we assume that we have permutation symmetry

if (iMethod == 2) then
  iOff = 1
  do iiSym=1,nSym
    iOrb = nRs1(iiSym)+nRs2(iiSym)+nRs3(iiSym)
    do jjSym=1,iiSym
      jOrb = nRs1(jjSym)+nRs2(jjSym)+nRs3(jjSym)

      if (Mul(iiSym,jjSym) == Dsym) then
        pINT1(iiSym) = iOff
        if (iiSym == jjSym) then
          iOff = iOff+nTri_Elem(iOrb)
        else
          iOff = iOff+iOrb*jOrb
        end if
      end if

    end do
  end do

  iOff = 1
  do iiSym=1,nSym
    iOrb = nRs1(iiSym)+nRs2(iiSym)+nRs3(iiSym)
    do jjSym=1,iiSym
      jOrb = nRs1(jjSym)+nRs2(jjSym)+nRs3(jjSym)

      ijSym = Mul(iiSym,jjSym)
      klSym = Mul(ijSym,DSym)

      ijNum = iTri(iiSym,jjSym)

      if (iiSym == jjSym) then
        ijOrb = nTri_Elem(iOrb)
      else
        ijOrb = iOrb*jOrb
      end if
      do kkSym=1,nsym
        kOrb = nRs1(kkSym)+nRs2(kkSym)+nRs3(kkSym)

        llSym = Mul(klSym,kkSym)
        lOrb = nRs1(llSym)+nRs2(llSym)+nRs3(llSym)

        if (llSym > kkSym) cycle

        klNum = iTri(kkSym,llSym)
        if (klNum > ijNum) cycle

        if (kkSym == llSym) then
          klOrb = nTri_Elem(kOrb)
        else
          klOrb = kOrb*lOrb
        end if
        ip = iiSym+nSym*((jjSym-1)+nSym*(kkSym-1))
        if (ijNum == klNum) then
          iPlus = nTri_Elem(ijOrb)
        else
          iPlus = ijOrb*klOrb
        end if

        if (iPlus > 0) pINT2(ip) = iOff

        iOff = iOff+iPlus

      end do
    end do
  end do
end if

call Get_iArray('nFro',nFro,nSym)
do i=1,nSym
  if (nFro(i) /= 0) then
    call WarningMessage(2,'MCLR module cannot handle frozen orbitals!')
    call Abend()
  end if
end do
call Put_iArray('nFroPT',nFro,nSym)
if (PT2) call Put_iArray('nIsh',nIsh,nSym)
call Get_iArray('nDel',nDel,nSym)
do i=1,nSym
  if (nDel(i) /= 0) then
    call WarningMessage(2,'MCLR module cannot handle deleted orbitals!')
    call Abend()
  end if
end do
call Put_iArray('nDelPT',nDel,nSym)

!call iWrtMa(pINT2,64,8,64,8)

end subroutine SetUp_MCLR
