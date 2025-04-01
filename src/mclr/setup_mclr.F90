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

use MCLR_Data, only: pInt1, pInt2
use MCLR_Data, only: nNA, n2Dens, nDens, nCMO, nDensC, nDens2, ipMatLT, ipMat, ipCM, ipMatBA, ipMO, nA, nB, n1Dens, nMBA
use input_mclr, only: nSym, TimeDep, iMethod, PT2, nAsh, nBas, nDel, nFro, nIsh, nOrb, nRS1, nRS2, nRs3

implicit none
integer dsym
! for the integrals needed in sigma gen
integer, external :: iPntSO
integer ip, nn, nBmx, iS, jS, lS, klS, kS, ijS, ipP, iExt0, iExt1, iExt2, iExt3, iInt4, iExt4, i1, iInt0, iInt1, iInt2, iInt3, &
        mATAB, iOff, iiSym, iOrb, jjSym, jOrb, ijSym, klSym, ijNum, ijOrb, kkSym, kOrb, llSym, lOrb, klNum, klOrb, iPlus, nDensLT
! Statement function
integer i, j, iTri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
ip = 1
call ICopy(8**3,[0],0,ipMO,1)

nn = 0
nbmx = 0
do iS=1,nsym
  nA(iS) = nn
  nB(is) = nAsh(is)+nIsh(is)
  nn = nn+nAsh(iS)
  nbmx = max(nbmx,nBas(iS))
end do

call Set_nbmx(nbmx)

nna = nn
n2Dens = itri(nnA**2,nnA**2)

n1Dens = nnA*nnA
ip = 1
do kS=1,nSym
  do lS=1,nsym
    klS = ieor(lS-1,kS-1)+1
    do jS=1,nSym
      do iS=1,nSym
        ijS = ieor(iS-1,jS-1)+1
        ipp = nOrb(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
        if ((ieor(ijS-1,klS-1)+1 == DSym) .and. (ipp /= 0)) then
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
call iCopy(8**2,[0],0,ipMat,1)
call iCopy(8**2,[0],0,ipMatLT,1)
call iCopy(8,[0],0,ipCM,1)

do jS=1,nSym
  do iS=1,js
    if (ieor(iS-1,jS-1) == DSym-1) then
      if (is < js) then
        ipMatLT(jS,iS) = ipntso(js-1,is-1,2**(dsym-1),nBas)+1
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
        nDensLT = nDensLT+nBas(iS)*(nBas(is)+1)/2
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
ndens2 = 0
matab = 1
do jS=1,nSym
  do iS=1,nSym
    if (ieor(iS-1,jS-1) == DSym-1) then
      ipMatba(is,js) = matab
      ipMat(jS,iS) = nDens2+1
      nDens2 = nDens2+nBas(iS)*nBas(jS)
      matab = matab+nash(js)*nOrb(iS)
    end if
  end do
end do
ndens = ndens2

! To begin with we assume that we have permutation symmetry

if (iMethod == 2) then
  iOff = 1
  do iiSym=1,nSym
    iOrb = nRs1(iiSym)+nRs2(iiSym)+nRs3(iiSym)
    do jjSym=1,iiSym
      jOrb = nRs1(jjSym)+nRs2(jjSym)+nRs3(jjSym)

      if (ieor(iiSym-1,jjSym-1)+1 == Dsym) then
        pINT1(iiSym) = iOff
        if (iiSym == jjSym) then
          iOff = iOff+iOrb*(iOrb+1)/2
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

      ijSym = ieor(iiSym-1,jjSym-1)+1
      klSym = ieor(ijSym-1,DSym-1)+1

      ijNum = iiSym*(iiSym+1)/2+jjSym

      if (iiSym == jjSym) then
        ijOrb = iOrb*(iOrb+1)/2
      else
        ijOrb = iOrb*jOrb
      end if
      do kkSym=1,nsym
        kOrb = nRs1(kkSym)+nRs2(kkSym)+nRs3(kkSym)

        llSym = ieor(klSym-1,kkSym-1)+1
        lOrb = nRs1(llSym)+nRs2(llSym)+nRs3(llSym)

        if (llSym > kkSym) cycle

        klNum = kkSym*(kkSym+1)/2+llSym
        if (klNum > ijNum) cycle

        if (kkSym == llSym) then
          klOrb = kOrb*(kOrb+1)/2
        else
          klOrb = kOrb*lOrb
        end if
        ip = iiSym+nSym*((jjSym-1)+nSym*(kkSym-1))
        if (ijNum == klNum) then
          iPlus = ijOrb*(ijOrb+1)/2
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
