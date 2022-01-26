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

use Basis_Info, only: dbsc, nBas, ncnttp
use Symmetry_Info, only: nIrrep
use DKH_Info, only: CLightAU
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
#include "rinfo.fh"
#include "print.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAa, iAaf, iAngr, iAux, iAux2, iBas, iBu, iBu2, iBu4, iBu6, iCmm1, iCmm2, icnt, iCnttp, iComp, iCr, idbg, &
                     iE, iEig, iEig4, iEigf, iEv2, iEv4, iEw, iEw4, iExp, i_f, if2, if2a, ifa, iG, iG2, iH, iH_nr, iH_temp, iip1, &
                     iK, iOpt, iOve, iP, ip1, ip2, ip4, ip5, ip6, ip7, ipaddr(3), iPrint, ipV, ipVf, ipVp, iRC, iRevt, iRevtf, &
                     iRout, iRr, iRrf, iScpV, iScVp, iSinv, iSinvf, iSize, iSizea, iSizeab, iSizeb, iSizec, iSizep, iSmlbl, iSS, &
                     iSyma, iSymb, iTt, iV, iVp, iVpf, ixyz, jExp, k, k1, k1a, k1b, k2, k2a, k2b, kAng, kC, kCof, kCofi, kCofj, &
                     kExp, kExpi, kExpj, kh, L, Len_, LenInt, LenIntf, LenIntf1, lOper, Lu_One, n, na, nb, ncomp, nSym
real(kind=wp) :: eps, rCofi, rCofj, rExpi, rExpj, rI, rNorm, Sum_, VELIT
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
if (iRC /= 0) Go To 9999

call OneBas('PRIM')
call Get_iArray('nBas_Prim',nbas,nSym)

if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Primitive basis fcns',(nBas(i),i=0,nSym-1)
end if

! Allocate memory for relativistic part

LenIntf = 0
LenIntf1 = 0
do L=0,nSym-1
  n = nBas(L)
  LenIntf = LenIntf+n*n+4
  LenIntf1 = LenIntf1+n+4
end do

call GetMem('Eigf    ','ALLO','REAL',iEigf,LenIntf)
call GetMem('Sinvf   ','ALLO','REAL',iSinvf,LenIntf)
call GetMem('Revtf   ','ALLO','REAL',iRevtf,LenIntf)
call GetMem('Aaf     ','ALLO','REAL',iAaf,LenIntf1)
call GetMem('Rrf     ','ALLO','REAL',iRrf,LenIntf1)
call dcopy_(LenIntf,[Zero],0,Work(iEigf),1)
call dcopy_(LenIntf,[Zero],0,Work(iSinvf),1)
call dcopy_(LenIntf,[Zero],0,Work(iRevtf),1)
call dcopy_(LenIntf1,[Zero],0,Work(iAaf),1)
call dcopy_(LenIntf1,[Zero],0,Work(iRrf),1)

VELIT = CLightAU
iSizep = 0
do L=0,nSym-1
  iSizep = iSizep+nBas(L)*(nBas(L)+1)/2
end do
if (iPrint >= 10) write(u6,*) ' iSizep',iSizep

call GetMem('Kin     ','ALLO','REAL',iK,iSizep+4)
call GetMem('SS      ','ALLO','REAL',iSS,iSizep+4)
call GetMem('V       ','ALLO','REAL',iV,iSizep+4)
call GetMem('pVp     ','ALLO','REAL',ipVp,iSizep+4)

if (iprint >= 20) write(u6,*) '  indices',iss,ik,iv,ipvp
Label = 'Mltpl  0'
iComp = 1
iOpt = 0
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
!write(u6,'(8f9.4)') (Work(iSS+k),k=0,iSizep-1)
nComp = 1
ipaddr(1) = iSS
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
Label = 'Attract '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'Kinetic '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iK),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'pVp     '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(ipVp),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

nComp = 3
call GetMem('lOper1  ','ALLO','INTE',ip2,nComp)
call GetMem('ipl     ','ALLO','INTE',ipV,nComp)
ixyz = 1
iWork(ip2) = 2**IrrFnc(ixyz)
ixyz = 2
iWork(ip2+1) = 2**IrrFnc(ixyz)
ixyz = 4
iWork(ip2+2) = 2**IrrFnc(ixyz)

call GetMem('iiVpf   ','ALLO','INTE',iVpf,nComp)
call GetMem('iipVf   ','ALLO','INTE',ipVf,nComp)

! Read pV matrix , elements <iSyma|pV|iSymb>, iSymb <= iSyma !

do iComp=1,nComp

  Label = 'pV      '
  iRC = -1
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)
  call GetMem('pV      ','ALLO','REAL',iWork(ipV+iComp-1),LenInt+4)
  ip4 = iWork(ipV+iComp-1)
  call RdOne(iRC,iOpt,Label,iComp,Work(ip4),iSmlbl)

end do

!call PrMtrx('pV      ',iWork(ip2),nComp,iWork(ipV),Work)

! Read Vp matrix , elements <iSyma|Vp|iSymb>, iSymb <= iSyma !

! Read Vp matrix

Label = 'Vp      '

call GetMem('ipk     ','ALLO','INTE',iVp,nComp)

do iComp=1,nComp

  iRC = -1
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)
  call GetMem('Vp      ','ALLO','REAL',iWork(iVp+iComp-1),LenInt+4)
  ip5 = iWork(iVp+iComp-1)
  call RdOne(iRC,iOpt,Label,iComp,Work(ip5),iSmlbl)
end do

!call PrMtrx(Label,iWork(ip2),nComp,iWork(iVp),Work)


! Build the whole matrix <iSyma|pV|iSymb> lower triangel
! and put into a Work(ip6+ip1+k) and the upper triange
! and put into a Work(ip6+LenInt+iip1+k).
! For 'Vp' matrix the same , p7 instead of p6

do iComp=1,nComp

  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)

  call GetMem('pVf     ','ALLO','REAL',iWork(ipVf+iComp-1),2*LenInt+4)
  call GetMem('Vpf     ','ALLO','REAL',iWork(iVpf+iComp-1),2*LenInt+4)

  ip4 = iWork(ipV+iComp-1)
  ip5 = iWork(iVp+iComp-1)
  ip6 = iWork(ipVf+iComp-1)
  ip7 = iWork(iVpf+iComp-1)
  ip1 = 0
  iip1 = 0

  do iSyma=0,nSym-1
    na = nBas(iSyma)
    if (na == 0) go to 81

    do iSymb=0,iSyma
      nb = nBas(iSymb)
      if (nb == 0) go to 82
      if (iand(iSmLbl,2**ieor(iSyma,iSymb)) == 0) Go to 82
      if (iSyma == iSymb) then
        iSize = na*(na+1)/2

        ! Elements <iSyma|pV|iSyma> and <iSyma|Vp|iSyma>

        do k=0,iSize-1
          Work(ip6+ip1+k) = Work(ip4+ip1+k)
          Work(ip7+ip1+k) = Work(ip5+ip1+k)
        end do
        ip1 = ip1+na*(na+1)/2
      else
        iSize = na*nb

        ! Elements <iSyma|pV|iSymb> and <iSymb|pV|iSyma> = -<iSyma|Vp|iSymb>
        ! AND
        ! Elements <iSyma|Vp|iSymb> and <iSymb|Vp|iSyma> = -<iSyma|pV|iSymb>

        do k=0,iSize-1
          Work(ip6+ip1+k) = Work(ip4+ip1+k)
          Work(ip7+ip1+k) = Work(ip5+ip1+k)
          Work(ip6+LenInt+4+iip1+k) = -Work(ip5+ip1+k)
          Work(ip7+LenInt+4+iip1+k) = -Work(ip4+ip1+k)
        end do
        ip1 = ip1+na*nb
        iip1 = iip1+na*nb
      end if

82    continue
    end do
81  continue
  end do

  !write(u6,*) 'Created  whole matrix (NAxNB) PV '
  !
  !write(u6,*) 'pV matrix <iSyma|pV|iSymb> '
  !write(u6,*)
  !write(u6,*) (Work(iWork(ipVf+iComp-1)+i),i=0,LenInt-1)
  !write(u6,*)
  !write(u6,*)
  !write(u6,*) 'pV matrix <iSymb|pV|iSyma> if iSyma /= iSymb otherwise zero matrix '
  !write(u6,*)
  !write(u6,*) (Work(iWork(ipVf+iComp-1)+LenInt+4+i),i=0,LenInt-1)
  !
  ! Create the whole matrix (NAxNB) VP
  !
  !write(u6,*) 'Create the whole matrix (NAxNB) VP '
  !
  !write(u6,*)
  !write(u6,*) 'Vp matrix <iSyma|Vp|iSymb> '
  !write(u6,*)
  !write(u6,*) (Work(iWork(iVpf+iComp-1)+i),i=0,LenInt-1)
  !write(u6,*)
  !write(u6,*)
  !write(u6,*)
  !write(u6,*) 'Vp matrix <iSymb|Vp|iSyma> if iSyma /= iSymb otherwise zero matrix '
  !write(u6,*)
  !write(u6,*)
  !write(u6,*) (Work(iWork(iVpf+iComp-1)+LenInt+4+i),i=0,LenInt-1)

  ! End of iComp loop

end do

! Main loop  1
eps = 1.0e-10_wp
k = 0
L = 0
k1 = 0
k2 = 0
do L=0,nSym-1

  n = nBas(L)
  iSize = n*(n+1)/2
  !AJS protection against zero dimension representation
  if (iSize == 0) goto 9
  !AJS put zeroes
  Len_ = 4*(iSize+4)+5*(n*n+4)+5*(n+4)
  call GetMem('Scratch ','ALLO','REAL',iCr,Len_+4)
  call dcopy_(Len_,[Zero],0,Work(iCr),1)
  call GetMem('Scratch ','FREE','REAL',iCr,Len_+4)

  ! Allocate

  call GetMem('Bu      ','ALLO','REAL',iBu,isize+4)
  call GetMem('P       ','ALLO','REAL',iP,isize+4)
  call GetMem('G       ','ALLO','REAL',iG,isize+4)
  call GetMem('Ev2     ','ALLO','REAL',iEv2,isize+4)
  call GetMem('Eig     ','ALLO','REAL',iEig,n*n+4)
  call GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
  call GetMem('Revt    ','ALLO','REAL',iRevt,n*n+4)
  call GetMem('Aux     ','ALLO','REAL',iAux,n*n+4)
  call GetMem('Ove     ','ALLO','REAL',iOve,n*n+4)
  call GetMem('Ew      ','ALLO','REAL',iEw,n+4)
  call GetMem('E       ','ALLO','REAL',iE,n+4)
  call GetMem('Aa      ','ALLO','REAL',iAa,n+4)
  call GetMem('Rr      ','ALLO','REAL',iRr,n+4)
  call GetMem('Tt      ','ALLO','REAL',iTt,n+4)

  ! Debug output on unit idbg
  if (IfTest) then
    idbg = 41
  else
    idbg = -1
  end if


  ! call to package relsewb
  call SCFCLI2(idbg,eps,Work(iSS+k),Work(iK+k),Work(iV+k),Work(ipVp+k),n,iSize,VELIT,Work(iBu),Work(iP),Work(iG),Work(iEv2), &
               Work(iEig),Work(iSinv),Work(iRevt),Work(iAux),Work(iOve),Work(iEw),Work(iE),Work(iAa),Work(iRr),Work(iTt))

  call dcopy_(n*n+4,Work(iEig),1,Work(iEigf+k1),1)
  call dcopy_(n*n+4,Work(iSinv),1,Work(iSinvf+k1),1)
  call dcopy_(n*n+4,Work(iRevt),1,Work(iRevtf+k1),1)
  call dcopy_(n+4,Work(iAa),1,Work(iAaf+k2),1)
  call dcopy_(n+4,Work(iRr),1,Work(iRrf+k2),1)

  call GetMem('Bu      ','FREE','REAL',iBu,isize+4)
  call GetMem('P       ','FREE','REAL',iP,isize+4)
  call GetMem('G       ','FREE','REAL',iG,isize+4)
  call GetMem('Ev2     ','FREE','REAL',iEv2,isize+4)
  call GetMem('Eig     ','FREE','REAL',iEig,n*n+4)
  call GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
  call GetMem('Revt    ','FREE','REAL',iRevt,n*n+4)
  call GetMem('Aux     ','FREE','REAL',iAux,n*n+4)
  call GetMem('Ove     ','FREE','REAL',iOve,n*n+4)
  call GetMem('Ew      ','FREE','REAL',iEw,n+4)
  call GetMem('E       ','FREE','REAL',iE,n+4)
  call GetMem('Aa      ','FREE','REAL',iAa,n+4)
  call GetMem('Rr      ','FREE','REAL',iRr,n+4)
  call GetMem('Tt      ','FREE','REAL',iTt,n+4)
9 k = k+isize
  k1 = k1+n*n+4
  k2 = k2+n+4
end do

! Main loop  2
nComp = 3

do iComp=1,nComp

  ip6 = iWork(ipVf+iComp-1)
  ip7 = iWork(iVpf+iComp-1)
  ip1 = 0
  iip1 = 0
  eps = 1.0e-10_wp
  kh = 0
  k1a = 0
  k2a = 0
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)

  do iSyma=0,nSym-1

    na = nBas(iSyma)
    iSizea = na*(na+1)/2
    !AJS protection against zero dimension representation
    if (iSizea <= 0) goto 19
    !AJS

    k1b = 0
    k2b = 0

    do iSymb=0,nSym-1

      nb = nBas(iSymb)
      iSizeb = nb*(nb+1)/2
      if (iSizeb <= 0) goto 29
      iSizeab = na*nb
      if (iand(iSmlbl,2**ieor(iSyma,iSymb)) == 0) Go to 29

      call GetMem('ipVa    ','ALLO','REAL',ifa,iSizea)
      call GetMem('iVpa    ','ALLO','REAL',if2a,iSizea)

      if (iSyma == iSymb) then

        call dcopy_(iSizea,[Zero],0,Work(ifa),1)
        call dcopy_(iSizea,[Zero],0,Work(if2a),1)
        do k=0,iSizea-1
          Work(ifa+k) = Work(ip6+ip1+k)
          Work(if2a+k) = Work(ip7+ip1+k)
        end do
        ip1 = ip1+iSizea
      end if

      call GetMem('ifpV    ','ALLO','REAL',i_f,iSizeab)
      call GetMem('ifVp    ','ALLO','REAL',if2,iSizeab)
      call GetMem('ScpV    ','ALLO','REAL',iScpV,iSizeab)
      call GetMem('ScVp    ','ALLO','REAL',iScVp,iSizeab)
      call dcopy_(iSizeab,[Zero],0,Work(i_f),1)
      call dcopy_(iSizeab,[Zero],0,Work(if2),1)
      call dcopy_(iSizeab,[Zero],0,Work(iScpV),1)
      call dcopy_(iSizeab,[Zero],0,Work(iScVp),1)

      if (iSyma > iSymb) then

        do k=0,iSizeab-1
          Work(i_f+k) = Work(ip6+ip1+k)
          Work(if2+k) = Work(ip7+ip1+k)
        end do
        ip1 = ip1+na*nb

      end if

      if (iSyma < iSymb) then

        do k=0,iSizeab-1
          Work(i_f+k) = Work(ip6+iip1+LenInt+4+k)
          Work(if2+k) = Work(ip7+iip1+LenInt+4+k)
        end do
        call DCOPY_(iSizeab,Work(i_f),1,Work(iScpV),1)
        call DCOPY_(iSizeab,Work(if2),1,Work(iScVp),1)

        iip1 = iip1+na*nb
      end if

      ! Allocate

      call GetMem('Bu2     ','ALLO','REAL',iBu2,na*nb+4)
      call GetMem('Bu4     ','ALLO','REAL',iBu4,na*nb+4)
      call GetMem('Bu6     ','ALLO','REAL',iBu6,na*na+4)
      call GetMem('Bu      ','ALLO','REAL',iBu,iSizea+4)
      call GetMem('P       ','ALLO','REAL',iP,iSizea+4)
      call GetMem('G2      ','ALLO','REAL',iG2,na*nb+4)
      call GetMem('Ev4     ','ALLO','REAL',iEv4,isizea+4)
      call GetMem('Eig4    ','ALLO','REAL',iEig4,na*na+4)
      !call GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
      !call GetMem('Revt    ','ALLO','REAL',iRevt,n*n+4)
      call GetMem('Aux2    ','ALLO','REAL',iAux2,na*nb+4)
      call GetMem('Cmm1    ','ALLO','REAL',iCmm1,na*nb+4)
      call GetMem('Cmm2    ','ALLO','REAL',iCmm2,na*nb+4)
      call GetMem('Ew4     ','ALLO','REAL',iEw4,na+4)
      !call GetMem('E       ','ALLO','REAL',iE,n+4)
      !call GetMem('Tt      ','ALLO','REAL',iTt,n+4)

      call dcopy_(na*nb+4,[Zero],0,Work(iBu2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iBu4),1)
      call dcopy_(na*na+4,[Zero],0,Work(iBu6),1)
      call dcopy_(na*na+4,[Zero],0,Work(iEig4),1)
      call dcopy_(na+4,[Zero],0,Work(iEw4),1)
      call dcopy_(iSizea+4,[Zero],0,Work(iEv4),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iG2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iAux2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iCmm1),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iCmm2),1)

      ! Calculate the matrix elements of the
      ! following commutator :
      ! Revta*AAA*<a|[bp,V]|b>store in CMM1(iSyma,iSymb) matrix,

      call VPBMBPV(idbg,eps,Work(iEigf+k1a),Work(iEigf+k1b),Work(iRevtf+k1a),Work(iSinvf+k1a),Work(iSinvf+k1b),na,nb,iSizea, &
                   iSizeb,VELIT,Work(iAaf+k2a),Work(iAaf+k2b),Work(iRrf+k2a),Work(iRrf+k2b),Work(i_f),Work(if2),iSyma,iSymb, &
                   Work(iBu2),Work(iG2),Work(iAux2),Work(iCmm1),Work(iBu4),Work(iCmm2),Work(ifa),Work(if2a),Work(iScpV),Work(iScVp))

      call dcopy_(na*nb+4,[Zero],0,Work(iBu2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iBu4),1)
      call dcopy_(na*na+4,[Zero],0,Work(iBu6),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iG2),1)
      call dcopy_(na*nb+4,[Zero],0,Work(iAux2),1)
      call dcopy_(na*na+4,[Zero],0,Work(iEig4),1)

      ! call to package relsewc, BSS up to the fourth order in alpha

      call SCFCLI4(idbg,eps,Work(iSS+kh),Work(iK+kh),Work(iRevtf+k1a),Work(iSinvf+k1a),na,nb,iSizea,iSizeb,VELIT,Work(iAaf+k2a), &
                   Work(iAaf+k2b),iSyma,iSymb,Work(iCmm1),Work(iCmm2),Work(iEv4),Work(iBu2),Work(iBu6),Work(iEig4),Work(iEw4), &
                   Work(iP))

      ! Free a space

      call GetMem('Bu2     ','FREE','REAL',iBu2,na*nb+4)
      call GetMem('Bu4     ','FREE','REAL',iBu4,na*nb+4)
      call GetMem('Bu6     ','FREE','REAL',iBu6,na*na+4)
      call GetMem('Bu      ','FREE','REAL',iBu,iSizea+4)
      call GetMem('P       ','FREE','REAL',iP,iSizea+4)
      call GetMem('G2      ','FREE','REAL',iG2,na*nb+4)
      call GetMem('Ev4     ','FREE','REAL',iEv4,isizea+4)
      call GetMem('Eig4    ','FREE','REAL',iEig4,na*na+4)
      !call GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
      !call GetMem('Revt    ','FREE','REAL',iRevt,n*n+4)
      call GetMem('Aux2    ','FREE','REAL',iAux2,na*nb+4)
      call GetMem('Cmm1    ','FREE','REAL',iCmm1,na*nb+4)
      call GetMem('Cmm2    ','FREE','REAL',iCmm2,na*nb+4)
      call GetMem('Ew4     ','FREE','REAL',iEw4,na+4)
      !call GetMem('E       ','FREE','REAL',iE,n+4)
      !call GetMem('Tt      ','FREE','REAL',iTt,n+4)

      call GetMem('ifpV    ','FREE','REAL',i_f,iSizeab)
      call GetMem('ifVp    ','FREE','REAL',if2,iSizeab)
      call GetMem('ScpV    ','FREE','REAL',iScpV,iSizeab)
      call GetMem('ScVp    ','FREE','REAL',iScVp,iSizeab)

      call GetMem('ipVa    ','FREE','REAL',ifa,iSizea)
      call GetMem('iVpa    ','FREE','REAL',if2a,iSizea)

29    continue
      k1b = k1b+nb*nb+4
      k2b = k2b+nb+4
    end do
19  kh = kh+iSizea
    k1a = k1a+na*na+4
    k2a = k2a+na+4
  end do
end do

do iComp=1,nComp
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)
  call GetMem('pVf     ','FREE','REAL',iWork(ipVf+iComp-1),2*LenInt+4)
  call GetMem('Vpf     ','FREE','REAL',iWork(iVpf+iComp-1),2*LenInt+4)
end do

do iComp=1,nComp
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)
  call GetMem('Vp      ','FREE','REAL',iWork(iVp+iComp-1),LenInt+4)
end do
call GetMem('ipk     ','FREE','INTE',iVp,nComp)

do iComp=1,nComp
  iSmlbl = iWork(ip2+iComp-1)
  LenInt = n2Tri(iSmlbl)
  call GetMem('pV      ','FREE','REAL',iWork(ipV+iComp-1),LenInt+4)
end do

call GetMem('iiVpf   ','FREE','INTE',iVpf,nComp)
call GetMem('iipVf   ','FREE','INTE',ipVf,nComp)

call GetMem('ipl     ','FREE','INTE',ipV,nComp)
call GetMem('lOper1  ','FREE','INTE',ip2,nComp)

call GetMem('pVp     ','FREE','REAL',ipVp,iSizep+4)

! open arrays in contracted basis

iSizec = 0
do L=1,nSym
  iSizec = iSizec+nrBas(L)*(nrBas(L)+1)/2
end do
call GetMem('H       ','ALLO','REAL',iH,iSizec+4)
call FZero(Work(iH),iSizec+4)
call GetMem('H_nr    ','ALLO','REAL',iH_nr,iSizec+4)
call FZero(Work(iH_nr),iSizec+4)
call GetMem('H_temp  ','ALLO','REAL',iH_temp,iSizec+4)
call FZero(Work(iH_temp),iSizec+4)

! compute stripped non-relativistic H

Label = 'Kinetic '
call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

Label = 'Attract '
call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

call DaXpY_(iSizep+4,One,Work(iSS),1,Work(iV),1)

call dcopy_(4,[Zero],0,Work(iH_temp+iSizec),1)
call repmat(idbg,Work(iV),Work(iH_temp),.true.)

call GetMem('V       ','FREE','REAL',iV,iSizep+4)
call GetMem('SS      ','FREE','REAL',iSS,iSizep+4)

! Close ONEREL and re-open ONEINT

iOpt = 0
iRC = -1
call ClsOne(iRC,iOpt)
if (iRC /= 0) Go To 9999
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
if (iRC /= 0) Go To 9999

! Transform to contracted basis

! The Hamiltonian is now in Kin.

call repmat(idbg,Work(iK),Work(iH),.true.)

call GetMem('Kin     ','FREE','REAL',iK,iSizep+4)

iOpt = 0
iRC = -1
Label = 'OneHam 0'
call RdOne(iRC,iOpt,Label,1,Work(iH_nr),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iOpt = 0
iRC = -1

! final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
! where (s) is stripped and with full charge

call DaXpY_(iSizec+4,-One,Work(iH_temp),1,Work(iH),1)
call DaXpY_(iSizec+4,One,Work(iH_nr),1,Work(iH),1)

call Get_iArray('nBas',nBas,nSym)
if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
end if
Label = 'OneHam 0'
lOper = 1
nComp = 1
ipaddr(1) = iH
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)

! Replace 1-el Hamiltonian on ONEINT

iRC = -1
call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
Label = 'OneHam  '
call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
if (iRC /= 0) then
  write(u6,*) 'BSSInt: Error writing to ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

call GetMem('OneHam  ','FREE','REAL',iH,iSizec+4)
call GetMem('H_nr    ','FREE','REAL',iH_nr,iSizec+4)
call GetMem('H_temp  ','FREE','REAL',iH_temp,iSizec+4)

call GetMem('Eigf    ','FREE','REAL',iEigf,LenIntf)
call GetMem('Sinvf   ','FREE','REAL',iSinvf,LenIntf)
call GetMem('Revtf   ','FREE','REAL',iRevtf,LenIntf)
call GetMem('Aaf     ','FREE','REAL',iAaf,LenIntf1)
call GetMem('Rrf     ','FREE','REAL',iRrf,LenIntf1)
return

9999 continue
write(u6,*) ' *** Error in subroutine BSSint ***'
write(u6,*) '     Abend in subroutine OpnOne'
call Abend()

end subroutine BSSint
