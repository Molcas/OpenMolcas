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

subroutine NEMO_Opt1()

use Basis_Info
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "rinfo.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
integer nBas_Prim(0:7), nBas_cont(0:7), lOper(3)
parameter(MxMltPl=10)
integer ipMP((MxMltPl+1)*(MxMltPl+2)*(MxMltPl+3)/6),iSm((MxMltPl+1)*(MxMltPl+2)*(MxMltPl+3)/6)
integer ip(3), iSml(3)
character*8 Label
real*8, dimension(:), allocatable :: P_Matrix, MP_Matrix
dimension Length(1), nint(1)

iRout = 77
iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
! Save basis set info from contracted run

if (iprint >= 10) write(6,*) ' In NEMO_Opt1',ncnttp
kCof = 0
kAng = 0
kExp = 0
kC = 0

! Normalize coefficients

do iCnttp=1,nCnttp

!-- Make a check that no cartesian d or higher have been used.
!   The reason for this restriction is found in tr_prm_cnt.
!   If that routine is generalized, then remove this check.

  lSh = 0
  kShStr = dbsc(iCnttp)%iVal
  kShEnd = dbsc(iCnttp)%iVal+dbsc(iCnttp)%nVal-1
  do kSh=kShStr,kShEnd
    if (.not. Shells(kSh)%Transf .and. lSh >= 2) then
      call WarningMessage(2,'   NEMO Error')
      write(6,*)
      write(6,*)
      write(6,*) 'Error! The NEMO keyword does not work with cartesian d-functions or higher.'
      write(6,*) 'Request spherical functions to proceed.'
      write(6,*)
      write(6,*)
      call Quit(_RC_INPUT_ERROR_)
    end if
    lSh = lSh+1
  end do

!-- End check.

  do icnt=1,dbsc(iCnttp)%nCntr
    kC = kC+1
    !do iAngr=0,nAngr(icnt)
    do iAngr=0,nAngr(kC)
      !rI = iAngr+1.0d0+half
      rI = dble(iAngr)+One+Half
      kAng = kAng+1
      do iBas=1,nBasisr(kAng)
        Sum = Zero
        kExpi = kExp
        kCofi = kCof
        do iExp=1,nPrimr(kAng)
          kExpi = kExpi+1
          kCofi = kCofi+1
          rExpi = rExp(kExpi)
          !write(6,'(a11,f20.8)') ' Exponents',rExpi
          rCofi = rCof(kCofi)
          kExpj = kExp
          kCofj = kCof
          do jExp=1,nPrimr(kAng)
            kExpj = kExpj+1
            kCofj = kCofj+1
            rExpj = rExp(kExpj)
            rCofj = rCof(kCofj)
            Sum = Sum+rCofi*rCofj*(Two*sqrt(rExpi*rExpj)/(rExpi+rExpj))**rI
          end do
        end do
        rNorm = One/sqrt(Sum)
        if (iprint >= 10) write(6,*) ' rNorm',kAng,rNorm
        do iExp=1,nPrimr(kAng)
          rCof(kCof+iExp) = rCof(kCof+iExp)*rNorm
          if (iprint >= 10) then
            write(6,'(a24,f20.6)') ' normalized coefficients',rCof(kCof+iExp)
          end if
        end do
        kCof = kCof+nPrimr(kAng)
      end do
      kExp = kExp+nPrimr(kAng)
    end do
  end do
end do

if (iPrint >= 10) then
  i = 0
  do L=1,nrSym
    write(6,*) ' Irreducible representation',L
    do ibas=1,nrBas(L)
      i = i+1
      write(6,'(20i4)') i,icent(i),lnang(i),lmag(i)
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Close ONEINT and re-open ONEREL

call iCopy(8,nBas,1,nBas_Cont,1)

nSym = nIrrep
iOpt = 0
call ClsOne(iRC,iOpt)
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
if (iRC /= 0) Go To 9999

call OneBas('PRIM')
call Get_iArray('nBas_Prim',nBas_Prim,nIrrep)

if (iPrint >= 10) then
  write(6,'(a,8i5)') ' Symmetries          ',nSym
  write(6,'(a,8i5)') ' Primitive basis fcns',(nBas_Prim(i),i=0,nSym-1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read P_Matrix from ONEREL

nComp = 3
nLength_Tot = 0
do iComp=1,nComp
  iOpt = 1
  iRC = -1
  call iRdOne(iRC,iOpt,'P_matrix',iComp,Length,iSmLbl)
  if (iRC /= 0) then
    call WarningMessage(2,'Error reading length of P-Matrix')
    write(6,*) 'iComp=',iComp
    call Abend()
  end if
  iSml(iComp) = iSmLbl
  lOper(iComp) = 1
  ip(iComp) = 1+nLength_Tot
  nLength_Tot = nLength_Tot+Length(1)+4
end do

call mma_allocate(P_Matrix,nLength_Tot,label='P_Matrix')
call FZero(P_Matrix,nLength_Tot)

do iComp=1,nComp
  iOpt = 0
  iRC = -1
  iSmLbl = iSml(iComp)
  ip(iComp) = ip(iComp)
  call RdOne(iRC,iOpt,'P_matrix',iComp,P_Matrix(ip(iComp)),iSmLbl)
  if (iRC /= 0) then
    call WarningMessage(2,'Error reading P-Matrix')
    write(6,*) 'iComp=',iComp
    call Abend()
  end if
end do
if (iPrint >= 10) call PrMtrx('P_matrix',lOper,nComp,ip,P_Matrix)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read multipole integrals from ONEREL

nip = 0
nInt_Tot = 0
do iMltPl=0,MxMltPl
  write(Label,'(a,i2)') 'MLTPL ',iMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp
    iRC = -1
    iOpt = 1
    nInt = 0
    call iRdOne(iRC,iOpt,Label,iComp,nInt,iSmLbl)
    if (iRC /= 0) then
      if (iComp /= 1) then
        call WarningMessage(2,' Error reading length!')
        write(6,*) ' Label=',Label,' Comp=',iComp
        call Abend()
      else
        Go To 100
      end if
    end if
    nip = nip+1
    iSm(nip) = iSmLbl
    ipMP(nip) = 1+nInt_Tot
    nInt_Tot = nInt_Tot+nint(1)+4
  end do
end do
100 continue

call mma_allocate(MP_Matrix,nInt_Tot,label='MP_Matrix')

iip = 0
do iMltPl=0,MxMltPl
  write(Label,'(a,i2)') 'MLTPL ',iMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp
    iip = iip+1
    if (iip > nip) Go To 200
    iRC = -1
    iOpt = 0
    ipMP(iip) = ipMP(iip)
    iSmLbl = iSm(iip)
    call RdOne(iRC,iOpt,Label,iComp,MP_Matrix(ipMP(iip)),iSmLbl)
    if (iRC /= 0) then
      call WarningMessage(2,' Error reading integrals!')
      write(6,*) ' Label=',Label,' Comp=',iComp
      call Abend()
    end if
  end do
end do
200 continue
!                                                                      *
!***********************************************************************
!                                                                      *
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

call OneBas('PRIM')
!                                                                      *
!***********************************************************************
!                                                                      *
! Put the transformation matrix on RUNFILE

idbg = 0
call tr_prm_cnt(idbg,nBas_Cont,nBas_Prim)
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the P-matrix

nComp = 3
iOpt = 0
do iComp=1,nComp
  iRC = -1
  call WrOne(iRC,iOpt,'P_matrix',iComp,P_Matrix(ip(iComp)),iSml(iComp))
  if (iRC /= 0) then
    call WarningMessage(2,'Error reading P-Matrix')
    write(6,*) 'iComp=',iComp
    call Abend()
  end if
end do
call mma_deallocate(P_Matrix)
!                                                                      *
!***********************************************************************
!                                                                      *
! Process multipole integrals

iip = 0
do iMltPl=0,MxMltPl
  write(Label,'(a,i2)') 'PLTPL ',iMltpl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  do iComp=1,nComp
    iip = iip+1
    if (iip > nip) Go To 300
    iRC = -1
    iOpt = 0
    call WrOne(iRC,iOpt,Label,iComp,MP_Matrix(ipMP(iip)),iSm(iip))
    if (iRC /= 0) then
      call WarningMessage(2,' Error writing integrals!')
      write(6,*) ' Label=',Label,' Comp=',iComp
      call Abend()
    end if
  end do
end do
300 continue
call mma_deallocate(MP_Matrix)
!                                                                      *
!***********************************************************************
!                                                                      *
! And now change it back!

call OneBas('CONT')
!                                                                      *
!***********************************************************************
!                                                                      *
return

9999 continue
call WarningMessage(2,' *** Error in subroutine NEMO_Opt1 ***;     Abend in subroutine OpnOne or ClsOne')
call Abend()

end subroutine NEMO_Opt1
