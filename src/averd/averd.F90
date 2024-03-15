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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************

subroutine Averd(ireturn)
!-- Compute average density and corresponding natural orbitals. Two
!   possibilities exists, either construct average from input orbitals,
!   or from density matrices.
!
!   Author: Anders Ohrn.

use Averd_global, only: Wset
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, iB, iB1, iB2, iC, icomp, iD, iDs, iDt, iDummy(7,8), iErr, ind, indB, indS, indt, iO, iopt, iPrint, iset, &
                     irc, iSym, isyml, itBas, j, k, kaunt, kaunter, lsmat, Luinp, LuOut, nB, nBS, nBT, nOrb, nSet, nSym, ntot, ntot2
real(kind=wp) :: Dum, Dummy(1), Sqroot, Thr, Thro, ThrOcc, Wsum
logical(kind=iwp) :: PrOcc, PrEne, DensityBased
character(len=72) :: Title
character(len=7) :: Fname
character(len=10) :: OLabel
character(len=40) :: Titorb
character(len=128) :: OrbFile
character(len=4000) :: BsLbl
character(len=3) :: PLab
integer(kind=iwp), external :: IsFreeUnit
integer(kind=iwp), allocatable :: nBas(:)
real(kind=wp), allocatable :: AUX(:), CMO(:), Dao(:), Dtemp(:), DTmp(:), Occ(:), OccNat(:), Occs(:), Orbs(:), OrtoD(:), OrtoDt(:), &
                              S(:), Si(:), Sp(:), Ss(:), St(:), Trani(:), Trans(:), Vecs(:), Zeros(:)
#include "Molcas.fh"

!-- Banner.

ireturn = 99

!-- Define defaults and initialize.

call Init_ave(Title,iPrint,PrOcc,PrEne,DensityBased,ThrOcc,Dummy(1),iDummy(1,1))

!-- Read input.

nSet = 0
call Get_Averd_input(Title,iPrint,nSet,DensityBased,ThrOcc)
if (.not. allocated(Wset)) call mma_allocate(Wset,nSet,label='Wset')

!-- Read some information from RUNFILE.

call Get_iScalar('nSym',nSym)
call mma_allocate(nBas,nSym,label='nBas')
call Get_iArray('nBas',nBas,nSym)
itBas = 0
do iSym=1,nSym
  itBas = itBas+nBas(isym)
end do
call Get_cArray('Unique Basis Names',BsLbl,LenIn8*itBas)

!-- Some dimensions.

lsmat = 0
ntot = 0
ntot2 = 0
do i=1,nSym
  lsmat = lsmat+(nBas(i)*(nBas(i)+1))/2
  ntot = ntot+nBas(i)
  ntot2 = ntot2+nBas(i)**2
end do

!-- Read AO-basis overlap matrix.

call mma_allocate(S,lsmat+4,label='Overlap')
OLabel = 'Mltpl  0'
irc = 0
iopt = ibset(ibset(0,sNoOri),sNoNuc)
icomp = 1
isyml = 1
call RdOne(irc,iopt,OLabel,icomp,S,isyml)
if (iprint >= 99) then
  ind = 1
  do iSym=1,nSym
    call TriPrt('Overlap Matrix',' ',S(ind),nBas(iSym))
    ind = ind+nBas(iSym)*(nBas(iSym)+1)/2
  end do
end if

!-- Normalize weights.

Wsum = sum(Wset(:))
Wset(:) = Wset(:)/Wsum

!-- Print some Bla Bla...

if (iprint >= 2) then
  call Print_Input(Title,nSym,nBas,wSet,nSet)
end if

!-- Do the dirty work. Different paths for orbital- and density-based
!   averageing.

call mma_allocate(Dao,ntot2,'Density')
Dao(:) = Zero
if (.not. DensityBased) then
  Luinp = 10
  call mma_allocate(CMO,ntot2,label='Orbitals')
  call mma_allocate(Occ,ntot,label='Occ')
  do iset=1,nSet
    Fname = 'NAT001'
    write(Fname(4:6),'(i3.3)') iset
    ! Read orbital coefficients and occupation numbers.
    call RdVec(Fname,Luinp,'CO',nSym,nBas,nBas,CMO,Occ,Dummy,iDummy,Titorb,0,iErr)
    iC = 1
    iO = 1
    iD = 1
    ! Up-date average density matrix.
    do isym=1,nSym
      kaunter = 0
      do i=1,nBas(iSym)
        do j=1,nBas(iSym)
          do k=1,nBas(iSym)
            Dao(iD+kaunter) = Dao(iD+kaunter)+Wset(iSet)*Occ(iO+k-1)*CMO(iC+i+(k-1)*nBas(iSym)-1)*CMO(iC+j+(k-1)*nBas(iSym)-1)
          end do
          kaunter = kaunter+1
        end do
      end do
      iC = iC+nBas(isym)**2
      iD = iD+nBas(isym)**2
      iO = iO+nBas(isym)
    end do
    ! Print print print.
    if (iPrint >= 5) then
      ThrO = 1.0e-5_wp
      call Primo(Titorb,PrOcc,PrEne,ThrO,Dummy(1),nSym,nBas,nBas,BsLbl,Dummy,Occ,CMO,-1)
    end if
  end do
  call mma_deallocate(CMO)
  call mma_deallocate(Occ)
else
  call mma_allocate(Dtemp,lsmat,label='DensityT')
  Dtemp(:) = Zero
  do iset=1,nSet
    Fname = 'RUN001'
    write(Fname(4:6),'(i3.3)') iset
    call NameRun(Fname)
    ! Collect density from runfile.
    call mma_allocate(DTmp,lsmat,Label='DTmp')
    call Get_dArray('D1ao',Dtmp,lsmat)
    call DaxPy_(lsmat,Wset(iset),Dtmp,1,Dtemp,1)
    call mma_deallocate(DTmp)
  end do
  ! Square the density matrix.
  iDt = 1
  iDs = 1
  do iSym=1,nSym
    nB = nBas(iSym)
    call Dsq(Dtemp(iDt),Dao(iDs),1,nB,nB)
    iDt = iDt+nB*(nB+1)/2
    iDs = iDs+nB**2
  end do
  call mma_deallocate(Dtemp)
end if

call mma_deallocate(Wset)

!-- With the average density in store, lets orthogonalize (canonical),
!   then diagonalize to get natural orbitals.

indT = 1
indS = 1
indB = 1
call mma_allocate(Orbs,ntot2,label='NatOrbAcc')
call mma_allocate(Occs,ntot,label='NatOccAcc')
do iSym=1,nSym
  nBT = nBas(iSym)*(nBas(iSym)+1)/2
  nBS = nBas(iSym)**2
  call mma_allocate(Vecs,nBS,label='EigV')
  call mma_allocate(St,nBT,label='St')
  call mma_allocate(Si,nBT,label='Si')
  call mma_allocate(Ss,nBS,label='Ss')
  call mma_allocate(Sp,nBS,label='Sp')
  call mma_allocate(AUX,nBS,label='AUX')
  call mma_allocate(Trans,nBS,label='TransS')
  call mma_allocate(Trani,nBS,label='TransSi')
  call mma_allocate(OrtoD,nBS,label='OrthoDensS')
  call mma_allocate(OrtoDt,nBT,label='OrthoDensT')
  call mma_allocate(OccNat,nBas(iSym),label='Occs')
  St(:) = Zero
  Si(:) = Zero
  kaunter = 1
  do iB1=1,nBas(iSym)
    do iB2=1,nBas(iSym)
      Vecs(kaunter) = Zero
      if (iB1 == iB2) Vecs(kaunter) = One
      kaunter = kaunter+1
    end do
  end do
  call Jacob(S(indT),Vecs,nBas(iSym),nBas(iSym))
  do i=1,nBas(iSym)
    Sqroot = sqrt(S(indT+i*(i+1)/2-1))
    St(i*(i+1)/2) = Sqroot
    Si(i*(i+1)/2) = One/Sqroot
  end do
  call Square(St,Ss,1,nBas(iSym),nBas(iSym))
  call Square(Si,Sp,1,nBas(iSym),nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Vecs,nBas(iSym),Ss,nBas(iSym),Zero,AUX,nBas(iSym))
  call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One,AUX,nBas(iSym),Vecs,nBas(iSym),Zero,Trans,nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Vecs,nBas(iSym),Sp,nBas(iSym),Zero,AUX,nBas(iSym))
  call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One,AUX,nBas(iSym),Vecs,nBas(iSym),Zero,Trani,nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Trans,nBas(iSym),Dao(indS),nBas(iSym),Zero,AUX,nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,AUX,nBas(iSym),Trans,nBas(iSym),Zero,OrtoD,nBas(iSym))
  kaunter = 1
  do iB1=1,nBas(iSym)
    do iB2=1,nBas(iSym)
      Vecs(kaunter) = Zero
      if (iB1 == iB2) Vecs(kaunter) = One
      kaunter = kaunter+1
    end do
  end do
  kaunter = 1
  do i=1,nBas(iSym)
    do j=1,i
      OrtoDt(kaunter) = OrtoD(i+(j-1)*nBas(iSym))
      kaunter = kaunter+1
    end do
  end do
  call Jacob(OrtoDt,Vecs,nBas(iSym),nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Trani,nBas(iSym),Vecs,nBas(iSym),Zero,AUX,nBas(iSym))
  kaunt = 1
  kaunter = 1
  do i=1,nBas(iSym)
    do j=1,i
      if (i == j) then
        OccNat(kaunt) = OrtoDt(kaunter)
        kaunt = kaunt+1
      end if
      kaunter = kaunter+1
    end do
  end do
  call SortEig(OccNat,AUX,nBas(iSym),nBas(iSym),-1,.false.)
  call Add_Info('AVERD_OCC',OccNat,5,5)
  if (iPrint >= 5) then
    write(Titorb,'(a)') 'All average orbitals in this irrep.'
    Thr = -One
    call Primo(Titorb,.true.,.false.,Thr,Dum,1,nBas(iSym),nBas(iSym),BsLbl,Dummy,OccNat,AUX,-1)
  end if
  call dcopy_(nBS,AUX,1,Orbs(indS),1)
  call dcopy_(nBas(iSym),OccNat,1,Occs(indB),1)
  call mma_deallocate(Vecs)
  call mma_deallocate(St)
  call mma_deallocate(Si)
  call mma_deallocate(Ss)
  call mma_deallocate(Sp)
  call mma_deallocate(AUX)
  call mma_deallocate(Trans)
  call mma_deallocate(Trani)
  call mma_deallocate(OrtoD)
  call mma_deallocate(OrtoDt)
  call mma_deallocate(OccNat)
  indT = indT+nBT
  indS = indS+nBS
  indB = indB+nBas(iSym)
end do
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,*) '---Average natural orbital generation completed!---'
write(u6,*)

!-- Write average orbitals to a file with the same format as
!   SCF-orbitals. To the outfile are orbital energies added just
!   to make the NEMO happy, the numbers are just bosh!

write(u6,*)
write(u6,*)
write(u6,*) 'Average orbitals put on AVEORB'
write(u6,*)
write(u6,*) 'NB: Dummy orbital energies added to AVEORB for compatability reasons.'
write(u6,*) '    They have no physical meaning.'
call mma_allocate(Zeros,ntot,label='Zeros')
Zeros(:) = Zero
LuOut = 65
LuOut = IsFreeUnit(LuOut)
Title = 'Average Orbitals'
OrbFile = 'AVEORB'
Plab = 'COE'
call WrVec(OrbFile,LuOut,Plab,nSym,nBas,nBas,Orbs,Occs,Zeros,iDummy,Title)

!-- Say something about orbital occupation.

write(u6,*)
write(u6,*)
write(u6,'(a)') ' |  Average orbital occupation.'
write(u6,'(a)') ' |-----------------------------'
write(u6,'(a,es18.8)') ' |    Threshold: ',ThrOcc
write(u6,*)
nOrb = 0
iO = 0
do iSym=1,nSym
  do iB=1,nBas(iSym)
    if (Occs(iO+iB) < ThrOcc) cycle
    nOrb = nOrb+1
  end do
  write(u6,'(a,i2,a,i4)') '      Symmetry:',iSym,'   Number of orbitals below threshold:',nOrb
  iO = iO+nBas(iSym)
end do
write(u6,*)
call mma_deallocate(nBas)
call mma_deallocate(Zeros)
call mma_deallocate(Orbs)
call mma_deallocate(Occs)
call mma_deallocate(Dao)
call mma_deallocate(S)

!-- Good Bye.

ireturn = 0

return

end subroutine Averd
