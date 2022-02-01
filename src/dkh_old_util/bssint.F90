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

subroutine BSSint()

use Basis_Info, only: dbsc, nBas, nCnttp
use Symmetry_Info, only: nIrrep
use DKH_Info, only: CLightAU
use Data_Structures, only: Alloc_Alloc1DArray, Alloc_Alloc2DArray, Alloc1DArray_Type, Alloc2DArray_Type, DSBA_Type, Allocate_DSBA, &
                           Deallocate_DSBA, Free_Alloc1DArray, Free_Alloc2DArray
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
#include "rinfo.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAngr, iBas, icnt, iCnttp, iComp, idbg, iExp, iip1, iOpt, ip1, iPrint, iRC, iRout, iSize, iSizea, iSizeab, &
                     iSizeb, iSizec, iSmlbl, iSyma, iSymb, ixyz, jExp, kAng, kC, kCof, kCofi, kCofj, kExp, kExpi, kExpj, L, &
                     LenInt, lOper, Lu_One, n, na, nb, nBasMax, ncomp, nSym
real(kind=wp) :: rCofi, rCofj, rExpi, rExpj, rI, rNorm, Sum_, VELIT
type(DSBA_Type) :: Eigf, Kin, pVp, Revtf, Sinvf, SS, V
integer(kind=iwp), allocatable :: lOper1(:)
real(kind=wp), allocatable :: Aa(:), Aux(:), Bu(:), Bu2(:), Bu4(:), Cmm1(:), Cmm2(:), E(:), Eig(:), Ev(:), Ew(:), G(:), H(:), &
                              H_nr(:), H_temp(:), ifpV(:), ifVp(:), ipVa(:), iVpa(:), Ove(:), P(:), Revt(:), Rr(:), ScpV(:), &
                              ScVp(:), Sinv(:), Tt(:)
type(Alloc1DArray_Type), allocatable :: Aaf(:), Rrf(:), pV(:), Vp(:)
type(Alloc2DArray_Type), allocatable :: pVf(:), Vpf(:)
character(len=8) :: Label
#ifdef _DEBUGPRINT_
#define _TEST_ .true.
#else
#define _TEST_ .false.
#endif
logical(kind=iwp), parameter :: IfTest = _TEST_
integer(kind=iwp), external :: IrrFnc, n2Tri

iRout = 77
iPrint = nPrint(iRout)

! Save basis set info from contracted run

if (iprint >= 10) write(u6,*) ' In dkint',ncnttp
kCof = 0
kAng = 0
kExp = 0
kC = 0

! Normalize coefficients

do iCnttp=1,nCnttp
  do icnt=1,dbsc(iCnttp)%nCntr
    kC = kC+1
    !do iAngr=0,nAngr(icnt)
    do iAngr=0,nAngr(kC)
      !rI = iAngr+OneHalf
      rI = real(iAngr,kind=wp)+OneHalf
      kAng = kAng+1
      do iBas=1,nBasisr(kAng)
        Sum_ = Zero
        kExpi = kExp
        kCofi = kCof
        do iExp=1,nPrimr(kAng)
          kExpi = kExpi+1
          kCofi = kCofi+1
          rExpi = rExp(kExpi)
          !write(u6,'(a11,f20.8)') ' Exponents',rExpi
          rCofi = rCof(kCofi)
          kExpj = kExp
          kCofj = kCof
          do jExp=1,nPrimr(kAng)
            kExpj = kExpj+1
            kCofj = kCofj+1
            rExpj = rExp(kExpj)
            rCofj = rCof(kCofj)
            Sum_ = Sum_+rCofi*rCofj*(Two*sqrt(rExpi*rExpj)/(rExpi+rExpj))**rI
          end do
        end do
        rNorm = One/sqrt(Sum_)
        if (iprint >= 10) write(u6,*) ' rNorm',kAng,rNorm
        do iExp=1,nPrimr(kAng)
          rCof(kCof+iExp) = rCof(kCof+iExp)*rNorm
          if (iprint >= 10) then
            write(u6,'(a24,f20.6)') ' normalized coefficients',rCof(kCof+iExp)
          end if
        end do
        kCof = kCof+nPrimr(kAng)
      end do
      kExp = kExp+nPrimr(kAng)
    end do
  end do
end do

i = 0
if (iPrint >= 10) then
  do L=1,nrSym
    write(u6,*) ' Irreducible representation',L
    do ibas=1,nrBas(L)
      i = i+1
      write(u6,'(20i4)') i,icent(i),lnang(i),lmag(i)
    end do
  end do
end if

! Close ONEINT and re-open ONEREL

nSym = nIrrep
iOpt = 0
call ClsOne(iRC,iOpt)
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
if (iRC /= 0) then
  write(u6,*) ' *** Error in subroutine BSSint ***'
  write(u6,*) '     Abend in subroutine OpnOne'
  call Abend()
end if

call OneBas('PRIM')
call Get_iArray('nBas_Prim',nBas,nSym)

if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Primitive basis fcns',(nBas(i),i=0,nSym-1)
end if

! Allocate memory for relativistic part

call Allocate_DSBA(Eigf,nBas,nBas,nSym,label='Eigf')
call Allocate_DSBA(Sinvf,nBas,nBas,nSym,label='Sinvf')
call Allocate_DSBA(Revtf,nBas,nBas,nSym,label='Revtf')
call Alloc_Alloc1DArray(Aaf,[1,nSym],label='Aaf')
call Alloc_Alloc1DArray(Rrf,[1,nSym],label='Rrf')
nBasMax = 0
do L=1,nSym
  n = nBas(L-1)
  nBasMax = max(nBasMax,n)
  call mma_allocate(Aaf(L)%A,n,label='Aaf_i')
  call mma_allocate(Rrf(L)%A,n,label='Rrf_i')
end do

VELIT = CLightAU

call Allocate_DSBA(Kin,nBas,nBas,nSym,aCase='TRI',label='Kin')
call Allocate_DSBA(SS,nBas,nBas,nSym,aCase='TRI',label='SS')
call Allocate_DSBA(V,nBas,nBas,nSym,aCase='TRI',label='V')
call Allocate_DSBA(pVp,nBas,nBas,nSym,aCase='TRI',label='pVp')

Label = 'Mltpl  0'
iComp = 1
iOpt = 0
iOpt = 6 ! Do not read origin or nuclear contribution
iRC = -1
call RdOne(iRC,iOpt,Label,1,SS%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
!write(u6,'(8f9.4)') SS%A0
if (iPrint >= 20) call PrMtrx(Label,[lOper],1,[1],SS%A0)
Label = 'Attract '
iRC = -1
call RdOne(iRC,iOpt,Label,1,V%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'Kinetic '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Kin%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'pVp     '
iRC = -1
call RdOne(iRC,iOpt,Label,1,pVp%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

nComp = 3
call mma_allocate(lOper1,nComp,label='lOper1')
do iComp=1,nComp
  ixyz = 2**(iComp-1)
  lOper1(iComp) = 2**IrrFnc(ixyz)
end do

! Read pV matrix , elements <iSyma|pV|iSymb>, iSymb <= iSyma !

call Alloc_Alloc1DArray(pV,[1,nComp],label='pV')

Label = 'pV      '
do iComp=1,nComp

  iRC = -1
  iSmlbl = lOper1(iComp)
  LenInt = n2Tri(iSmlbl)
  call mma_allocate(pV(iComp)%A,LenInt,label='pV_i')
  call RdOne(iRC,iOpt,Label,iComp,pV(iComp)%A,iSmlbl)

end do

! Read Vp matrix , elements <iSyma|Vp|iSymb>, iSymb <= iSyma !

! Read Vp matrix

call Alloc_Alloc1DArray(Vp,[1,nComp],label='Vp')

Label = 'Vp      '
do iComp=1,nComp

  iRC = -1
  iSmlbl = lOper1(iComp)
  LenInt = n2Tri(iSmlbl)
  call mma_allocate(Vp(iComp)%A,LenInt,label='Vp_i')
  call RdOne(iRC,iOpt,Label,iComp,Vp(iComp)%A,iSmlbl)
end do

! Build the whole matrix <iSyma|pV|iSymb> lower triangle
! and put into pVf(iComp)%A(:,1) and the upper triangle
! and put into pVf(iComp)%A(:,2).
! For 'Vp' matrix the same, Vpf instead of pVf

call Alloc_Alloc2DArray(pVf,[1,nComp],label='pVf')
call Alloc_Alloc2DArray(Vpf,[1,nComp],label='Vpf')

do iComp=1,nComp

  iSmlbl = lOper1(iComp)
  LenInt = n2Tri(iSmlbl)

  call mma_allocate(pVf(iComp)%A,LenInt,2,label='pVf_i')
  call mma_allocate(Vpf(iComp)%A,LenInt,2,label='Vpf_i')

  ip1 = 0
  iip1 = 0

  do iSyma=1,nSym
    na = nBas(iSyma-1)
    if (na == 0) cycle

    do iSymb=1,iSyma
      nb = nBas(iSymb-1)
      if (nb == 0) cycle
      if (iand(iSmLbl,2**ieor(iSyma-1,iSymb-1)) == 0) cycle
      if (iSyma == iSymb) then
        iSize = na*(na+1)/2

        ! Elements <iSyma|pV|iSyma> and <iSyma|Vp|iSyma>

        pVf(iComp)%A(ip1+1:ip1+iSize,1) = pV(iComp)%A(ip1+1:ip1+iSize)
        Vpf(iComp)%A(ip1+1:ip1+iSize,1) = Vp(iComp)%A(ip1+1:ip1+iSize)
        ip1 = ip1+na*(na+1)/2
      else
        iSize = na*nb

        ! Elements <iSyma|pV|iSymb> and <iSymb|pV|iSyma> = -<iSyma|Vp|iSymb>
        ! AND
        ! Elements <iSyma|Vp|iSymb> and <iSymb|Vp|iSyma> = -<iSyma|pV|iSymb>

        pVf(iComp)%A(ip1+1:ip1+iSize,1) = pV(iComp)%A(ip1+1:ip1+iSize)
        Vpf(iComp)%A(ip1+1:ip1+iSize,1) = Vp(iComp)%A(ip1+1:ip1+iSize)
        pVf(iComp)%A(iip1+1:iip1+iSize,2) = -Vp(iComp)%A(ip1+1:ip1+iSize)
        Vpf(iComp)%A(iip1+1:iip1+iSize,2) = -pV(iComp)%A(ip1+1:ip1+iSize)
        ip1 = ip1+na*nb
        iip1 = iip1+na*nb
      end if

    end do
  end do

  !write(u6,*) 'Created  whole matrix (NAxNB) PV '
  !
  !write(u6,*) 'pV matrix <iSyma|pV|iSymb> '
  !write(u6,*)
  !write(u6,*) pVf(iComp)%A(:,1)
  !write(u6,*)
  !write(u6,*)
  !write(u6,*) 'pV matrix <iSymb|pV|iSyma> if iSyma /= iSymb otherwise zero matrix '
  !write(u6,*)
  !write(u6,*) pVf(iComp)%A(:,2)
  !
  ! Create the whole matrix (NAxNB) VP
  !
  !write(u6,*) 'Create the whole matrix (NAxNB) VP '
  !
  !write(u6,*)
  !write(u6,*) 'Vp matrix <iSyma|Vp|iSymb> '
  !write(u6,*)
  !write(u6,*) Vpf(iComp)%A(:,1)
  !write(u6,*)
  !write(u6,*)
  !write(u6,*) 'Vp matrix <iSymb|Vp|iSyma> if iSyma /= iSymb otherwise zero matrix '
  !write(u6,*)
  !write(u6,*) Vpf(iComp)%A(:,2)

  ! End of iComp loop

end do

call mma_allocate(Bu,nBasMax*(nBasMax+1)/2,label='Bu')
call mma_allocate(P,nBasMax*(nBasMax+1)/2,label='P')
call mma_allocate(Ev,nBasMax*(nBasMax+1)/2,label='Ev')
call mma_allocate(G,nBasMax**2,label='G')
call mma_allocate(Eig,nBasMax**2,label='Eig')
call mma_allocate(Sinv,nBasMax**2,label='Sinv')
call mma_allocate(Revt,nBasMax**2,label='Revt')
call mma_allocate(Aux,nBasMax**2,label='Aux')
call mma_allocate(Ove,nBasMax**2,label='Ove')
call mma_allocate(Ew,nBasMax,label='Ew')
call mma_allocate(E,nBasMax,label='E')
call mma_allocate(Aa,nBasMax,label='Aa')
call mma_allocate(Rr,nBasMax,label='Rr')
call mma_allocate(Tt,nBasMax,label='Tt')

! Main loop  1
do L=1,nSym

  n = nBas(L-1)
  iSize = n*(n+1)/2
  !AJS protection against zero dimension representation
  if (iSize == 0) cycle

  ! Debug output on unit idbg
  if (IfTest) then
    idbg = 41
  else
    idbg = -1
  end if

  ! call to package relsewb
  call SCFCLI2(idbg,SS%SB(L)%A1,Kin%SB(L)%A1,V%SB(L)%A1,pVp%SB(L)%A1,n,iSize,VELIT,Bu,P,G,Ev,Eig,Sinv,Revt,Aux,Ove,Ew,E,Aa,Rr,Tt)

  Eigf%SB(L)%A1(:) = Eig(1:n*n)
  Sinvf%SB(L)%A1(:) = Sinv(1:n*n)
  Revtf%SB(L)%A1(:) = Revt(1:n*n)
  Aaf(L)%A(:) = Aa(1:n)
  Rrf(L)%A(:) = Rr(1:n)

end do

call mma_deallocate(Bu)
call mma_deallocate(Sinv)
call mma_deallocate(Revt)
call mma_deallocate(Ove)
call mma_deallocate(E)
call mma_deallocate(Aa)
call mma_deallocate(Rr)
call mma_deallocate(Tt)

call mma_allocate(ipVa,nBasMax*(nBasMax+1)/2,label='ipVa')
call mma_allocate(iVpa,nBasMax*(nBasMax+1)/2,label='iVpa')
call mma_allocate(ifpV,nBasMax**2,label='ifpV')
call mma_allocate(ifVp,nBasMax**2,label='ifVp')
call mma_allocate(ScpV,nBasMax**2,label='ScpV')
call mma_allocate(ScVp,nBasMax**2,label='ScVp')
call mma_allocate(Bu2,nBasMax**2,label='Bu2')
call mma_allocate(Bu4,nBasMax**2,label='Bu4')
call mma_allocate(Cmm1,nBasMax**2,label='Cmm1')
call mma_allocate(Cmm2,nBasMax**2,label='Cmm2')

! Main loop  2
do iComp=1,nComp

  ip1 = 0
  iip1 = 0
  iSmlbl = lOper1(iComp)
  LenInt = n2Tri(iSmlbl)

  do iSyma=1,nSym

    na = nBas(iSyma-1)
    iSizea = na*(na+1)/2
    !AJS protection against zero dimension representation
    if (iSizea <= 0) cycle
    !AJS

    do iSymb=1,nSym

      nb = nBas(iSymb-1)
      iSizeb = nb*(nb+1)/2
      if (iSizeb <= 0) cycle
      iSizeab = na*nb
      if (iand(iSmlbl,2**ieor(iSyma-1,iSymb-1)) == 0) cycle

      if (iSyma == iSymb) then

        ipVa(1:iSizea) = pVf(iComp)%A(ip1+1:ip1+iSizea,1)
        iVpa(1:iSizea) = Vpf(iComp)%A(ip1+1:ip1+iSizea,1)
        ip1 = ip1+iSizea
      end if

      if (iSyma > iSymb) then

        ifpV(1:iSizeab) = pVf(iComp)%A(ip1+1:ip1+iSizeab,1)
        ifVp(1:iSizeab) = Vpf(iComp)%A(ip1+1:ip1+iSizeab,1)
        ip1 = ip1+na*nb

      else if (iSyma < iSymb) then

        ifpV(1:iSizeab) = pVf(iComp)%A(iip1+1:iip1+iSizeab,2)
        ifVp(1:iSizeab) = Vpf(iComp)%A(iip1+1:iip1+iSizeab,2)
        ScpV(1:iSizeab) = ifpV(1:iSizeab)
        ScVp(1:iSizeab) = ifVp(1:iSizeab)

        iip1 = iip1+na*nb

      else

        ifpV(1:iSizeab) = Zero
        ifVp(1:iSizeab) = Zero
        ScpV(1:iSizeab) = Zero
        ScVp(1:iSizeab) = Zero

      end if

      Bu2(1:na*nb) = Zero
      Bu4(1:na*nb) = Zero
      Ew(1:na) = Zero
      Ev(1:iSizea) = Zero
      G(1:na*nb) = Zero
      Aux(1:na*nb) = Zero
      Cmm1(1:na*nb) = Zero
      Cmm2(1:na*nb) = Zero

      ! Calculate the matrix elements of the
      ! following commutator :
      ! Revta*AAA*<a|[bp,V]|b>store in CMM1(iSyma,iSymb) matrix,

      call VPBMBPV(Eigf%SB(iSyma)%A2,Eigf%SB(iSymb)%A2,Revtf%SB(iSyma)%A2,Sinvf%SB(iSyma)%A2,Sinvf%SB(iSymb)%A2,na,nb,iSizea, &
                   Aaf(iSyma)%A,Aaf(iSymb)%A,Rrf(iSyma)%A,Rrf(iSymb)%A,ifpV,ifVp,iSyma,iSymb,Bu2,G,Aux,Cmm1,Bu4,Cmm2,ipVa,iVpa, &
                   ScpV,ScVp)

      Bu2(1:na*na) = Zero
      Aux(1:na*nb) = Zero
      Eig(1:na*na) = Zero

      ! call to package relsewc, BSS up to the fourth order in alpha

      call SCFCLI4(idbg,SS%SB(iSyma)%A1,Kin%SB(iSyma)%A1,Sinvf%SB(iSyma)%A2,na,nb,iSizea,VELIT,Cmm1,Cmm2,Ev,Bu2,Eig,Ew,P)

    end do
  end do
end do

call mma_deallocate(P)
call mma_deallocate(Ev)
call mma_deallocate(G)
call mma_deallocate(Eig)
call mma_deallocate(Aux)
call mma_deallocate(Ew)
call mma_deallocate(ipVa)
call mma_deallocate(iVpa)
call mma_deallocate(ifpV)
call mma_deallocate(ifVp)
call mma_deallocate(ScpV)
call mma_deallocate(ScVp)
call mma_deallocate(Bu2)
call mma_deallocate(Bu4)
call mma_deallocate(Cmm1)
call mma_deallocate(Cmm2)

call mma_deallocate(lOper1)

call Free_Alloc2DArray(pVf)
call Free_Alloc2DArray(Vpf)
call Free_Alloc1DArray(pV)
call Free_Alloc1DArray(Vp)
call Deallocate_DSBA(pVp)

! open arrays in contracted basis

iSizec = 0
do L=1,nSym
  iSizec = iSizec+nrBas(L)*(nrBas(L)+1)/2
end do
call mma_allocate(H,iSizec+4)
call mma_allocate(H_nr,iSizec+4)
call mma_allocate(H_temp,iSizec+4)
H(:) = Zero
H_nr(:) = Zero
H_temp(:) = Zero

! compute stripped non-relativistic H

Label = 'Kinetic '
call RdOne(iRC,iOpt,Label,1,SS%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

Label = 'Attract '
call RdOne(iRC,iOpt,Label,1,V%A0,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

V%A0(:) = V%A0+SS%A0

call repmat(idbg,V%A0,H_temp,.true.)

call Deallocate_DSBA(SS)
call Deallocate_DSBA(V)

! Close ONEREL and re-open ONEINT

iOpt = 0
iRC = -1
call ClsOne(iRC,iOpt)
if (iRC /= 0) then
  write(u6,*) ' *** Error in subroutine BSSint ***'
  write(u6,*) '     Abend in subroutine OpnOne'
  call Abend()
end if
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
if (iRC /= 0) then
  write(u6,*) ' *** Error in subroutine BSSint ***'
  write(u6,*) '     Abend in subroutine OpnOne'
  call Abend()
end if

! Transform to contracted basis

! The Hamiltonian is now in Kin.

call repmat(idbg,Kin%A0,H,.true.)

call Deallocate_DSBA(Kin)

iOpt = 0
iRC = -1
Label = 'OneHam 0'
call RdOne(iRC,iOpt,Label,1,H_nr,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iRC = -1

! final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
! where (s) is stripped and with full charge

H(:) = H(:)-H_temp+H_nr

call Get_iArray('nBas',nBas,nSym)
if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
end if
Label = 'OneHam 0'
lOper = 1
if (iPrint >= 20) call PrMtrx(Label,[lOper],1,[1],H)

! Replace 1-el Hamiltonian on ONEINT

iRC = -1
call WrOne(iRC,iOpt,Label,1,H,lOper)
Label = 'OneHam  '
call WrOne(iRC,iOpt,Label,1,H,lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error writing to ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

call mma_deallocate(H)
call mma_deallocate(H_nr)
call mma_deallocate(H_temp)

call Deallocate_DSBA(Eigf)
call Deallocate_DSBA(Sinvf)
call Deallocate_DSBA(Revtf)
call Free_Alloc1DArray(Aaf)
call Free_Alloc1DArray(Rrf)

return

contains

end subroutine BSSint
