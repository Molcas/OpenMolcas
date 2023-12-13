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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine GiveMeInfo(nntyp,natyp,BasCoo,iCon,nPrim,nBA,nCBoA,nBonA,Expo,Cont,nSh,nfSh,nSize,iPrint,nAtoms,MxAngqNr,Acc,nBas)

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Real_Spherical, only: ipSph, RSph
use Index_Functions, only: nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nAtoms, MxAngqNr
integer(kind=iwp), intent(out) :: nntyp, natyp(nAtoms), nBA(nAtoms), nCBoA(nAtoms,MxAngqNr), nBonA(nAtoms), nSh(nAtoms), &
                                  nfSh(nAtoms,MxAngqNr), nSize, nBas
integer(kind=iwp), allocatable, intent(out) :: iCon(:,:), nPrim(:)
real(kind=wp), allocatable, intent(out) :: Expo(:,:), Cont(:,:), BasCoo(:,:), Acc(:)
integer(kind=iwp) :: i, iAng, iAngSav, iBas, iCnt, iCnttp, iCount, iHowMuch, ii, ind, ind1, ind2, ind3, ioio, iPrim, iTemp, j, jj, &
                     jSum, k, kaunt, kaunta, kaunter, kaunterPrev, kauntSav, kk, kkk, krekna, krekna2, l, ll, MaxAng, MxPrCon, na, &
                     nACCSize, ndc, nDiff, nnaa, nshj, nSumma, nVarv
real(kind=wp), allocatable :: TEMP1(:), TEMP2(:)
logical(kind=iwp) :: DoRys
#include "warnings.h"

!----------------------------------------------------------------------*
! Initialize in order to read properly from the info file.             *
!----------------------------------------------------------------------*
call Seward_Init()

!----------------------------------------------------------------------*
! GetInf reads everything in the runfile and puts it in variables      *
! in modules.                                                          *
!----------------------------------------------------------------------*
nDiff = 0
DoRys = .false.
call GetInf(DoRys,nDiff)

!----------------------------------------------------------------------*
! Set nntyp.                                                           *
!----------------------------------------------------------------------*
nntyp = nCnttp

!----------------------------------------------------------------------*
! Compute what we came here for. iBasAng will contain nBas elements    *
! with integers, such that 1=s-orbitals, 2=p-orbitals, 3=d-orbitals... *
!----------------------------------------------------------------------*
!ii = 0 !ii is number of basis sets.
!do
!  ii = ii+1
!  if (dbsc(ii)%nCntr == 0) exit
!end do
!ii = ii-1
!if (ii == 0) then
!  write(u6,*)
!  write(u6,*) 'ERROR in GiveMeInfo. No atoms?'
!end if
ii = nntyp

kaunta = 0
kaunt = 0
kaunter = 0
krekna = 0
MaxAng = 0
do i=1,ii
  kauntSav = kaunt
  do ioio=1,dbsc(i)%nCntr
    krekna = krekna+1
    krekna2 = 0
    kaunt = kauntSav
    kaunterPrev = kaunter
    nBA(krekna) = dbsc(i)%nShells
    if (nBA(krekna) > MaxAng) MaxAng = nBA(krekna)
    do j=1,dbsc(i)%nShells
      kaunt = kaunt+1
      krekna2 = krekna2+1
      nCBoA(krekna,krekna2) = Shells(kaunt)%nBasis
      kaunter = kaunter+Shells(kaunt)%nBasis
    end do
    kaunta = kaunta+1
    nBonA(kaunta) = kaunter-kaunterPrev  !Number of bases on each atom used below.
  end do
end do
nBas = kaunter
call mma_allocate(BasCoo,3,nBas,label='BasCoo')

!----------------------------------------------------------------------*
! And now coordinates of each basis.                                   *
!----------------------------------------------------------------------*
kaunter = 0
kaunt = 0
do i=1,ii
  do j=1,dbsc(i)%nCntr
    kaunt = kaunt+1
    do kk=1,nBonA(kaunt)
      kaunter = kaunter+1
      BasCoo(:,kaunter) = dbsc(i)%Coor(:,j)
    end do
  end do
end do

!----------------------------------------------------------------------*
! Now get info regarding the contraction. Icon is an array that for    *
! each basis type contain n1+n2+...+nx elements where n1 is the number *
! of contracted basis functions of s-type, n2 the same number for      *
! p-type etc. The value of the first n1 elements is the number of      *
! primitive basis functions of s-type, etc. So a contraction 7s3p.4s1p *
! generates the vector 7,7,7,7,3. We also compute natyp and also       *
! collect all exponents and contraction coefficients. These are stored *
! dynamically and then we return the pointers only.                    *
!----------------------------------------------------------------------*
MxPrCon = 0
kaunt = 0
do i=1,ii
  kaunter = 0
  do k=1,dbsc(i)%nShells
    kaunt = kaunt+1
    kaunter = kaunter+Shells(kaunt)%nBasis
  end do
  MxPrCon = max(MxPrCon,kaunter)
end do

call mma_allocate(Icon,nAtoms,MxPrCon,label='Icon')
kaunt = 0
do i=1,ii
  kaunter = 0
  do k=1,dbsc(i)%nShells
    kaunt = kaunt+1
    do ll=1,Shells(kaunt)%nBasis
      kaunter = kaunter+1
      Icon(i,kaunter) = Shells(kaunt)%nExp
    end do
  end do
end do

ndc = 0
iAngSav = 1
nSize = 0
kaunt = 0
do kk=1,ii !Just to get size of vector
  do kkk=1,dbsc(kk)%nShells
    kaunt = kaunt+1
    nSize = nSize+Shells(kaunt)%nBasis*Shells(kaunt)%nExp
  end do
end do
call mma_allocate(Expo,nntyp,nSize,label='Exponents')
call mma_allocate(Cont,nntyp,nSize,label='ContrCoef')
Expo(:,:) = Zero
Cont(:,:) = Zero

do iCnttp=1,nCnttp  !Here we set NaTyp.
  jSum = 0
  iTemp = 0
  nVarv = dbsc(iCnttp)%nShells
  nSh(iCnttp) = nVarv
  do iCnt=1,dbsc(iCnttp)%nCntr
    ndc = ndc+1
    iTemp = iTemp+dc(ndc)%nStab
  end do
  NaTyp(iCnttp) = iTemp
  do iAng=0,nVarv-1  !And in this loop we get hold of the contraction coefficients and the exponents.
    iCount = iAng+iAngSav
    iPrim = Shells(iCount)%nExp
    iBas = Shells(iCount)%nBasis
#   ifdef _DEBUGPRINT_
    call RecPrt('Exp',' ',Shells(iCount)%Exp,iPrim,1)
    call RecPrt('Cff',' ',Shells(iCount)%pCff,iPrim,iBas)
#   endif
    nfSh(iCnttp,iAng+1) = iBas
    do i=1,iBas
      Expo(iCnttp,jSum+1:jSum+iPrim) = Shells(iCount)%Exp(1:iPrim)
      Cont(iCnttp,jSum+1:jSum+iPrim) = Shells(iCount)%pCff(1:iPrim,i)
      jSum = jSum+iPrim
    end do
  end do
  iAngSav = iAngSav+iAng
end do
if (iPrint >= 30) then
  write(u6,*) 'Exp.'
  write(u6,'(10G13.4)') Expo(:,:)
  write(u6,*) 'Contr.'
  write(u6,'(10G13.4)') Cont(:,:)
end if

!----------------------------------------------------------------------*
! Construct the nPrim vector.                                          *
!----------------------------------------------------------------------*
call mma_allocate(nPrim,nbas,label='nPrim')
iBas = 0
do i=1,nntyp
  na = natyp(i)
  do j=1,na
    ind = 0
    nshj = nsh(i)
    do k=1,nshj
      nnaa = nfsh(i,k)
      do l=1,nnaa
        iBas = iBas+1
        ind = ind+1
        nPrim(iBas) = iCon(i,ind)
      end do
    end do
  end do
end do

! Then since overlap integrations are in cartesian coordinates while
! the AO-basis is spherical, we need transformation matrix for this.
! To our great joy, old reliable Seward computes this matrix of any
! order (within Molcas limits). Due to conflicting order conventions,
! some numbers gymnastics are required.

MaxAng = MaxAng-1
nSize = (2*MaxAng+1)*nTri_Elem1(MaxAng)
nACCSize = 0
do i=2,MaxAng
  nACCSize = nACCSize+(2*i+1)*nTri_Elem1(i)
end do
call mma_allocate(TEMP1,nSize,Label='TEMP1')
call mma_allocate(TEMP2,nSize,Label='TEMP2')
call mma_allocate(Acc,nACCSize,label='AccTransa')

nSumma = 0
do i=2,MaxAng
  ind1 = nTri_Elem1(i)
  ind2 = 2*i+1
  iHowMuch = ind1*ind2
  TEMP1(1:iHowMuch) = RSph(ipSph(i):ipSph(i)+iHowMuch-1)
  ind3 = 1
  do jj=1,ind1
    call dcopy_(ind2,TEMP1(jj),ind1,TEMP2(ind3),1)
    ind3 = ind3+ind2
  end do
  !call recprt('FFF',' ',TEMP1,nTri_Elem1(i),2*i+1)
  !call recprt('GGG',' ',TEMP2,ind2,ind1)
  Acc(nSumma+1:nSumma+iHowMuch) = Temp2(1:iHowMuch)
  nSumma = nSumma+iHowMuch
end do

call mma_deallocate(TEMP1)
call mma_deallocate(TEMP2)
!----------------------------------------------------------------------*
! Make deallocations. They are necessary because of the getinf.        *
!----------------------------------------------------------------------*
call ClsSew()

return

end subroutine GiveMeInfo
