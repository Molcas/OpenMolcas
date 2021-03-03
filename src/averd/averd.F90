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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "mxdm.fh"
#include "mxave.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAUX, iB, iB1, iB2, iC, iCMO, icomp, iD, iDao, iDs, iDt, iDtemp, iDummy(7,8), iErr, ind, indB, indS, indt, &
                     iO, iOcc, iOccNat, iOccs, iopt, iOrbs, iOrtoD, iOrtoDt, iPrint, iS, iset, iSi, irc, iSp, iSs, iSt, iSym, &
                     isyml, itBas, iTrani, iTrans, iVecs, iZero, j, k, kaunt, kaunter, lsmat, Luinp, LuOut, nB, nBas(MxSym), nBS, &
                     nBT, nOrb, Nset, nSym, ntot, ntot2
real(kind=wp) :: Dum, Dummy(1), Sqroot, Thr, Thro, ThrOcc, Wset(MxSets), Wsum
logical(kind=iwp) :: PrOcc, PrEne, DensityBased
character(len=72) :: Title
character(len=7) :: Fname
character(len=10) :: OLabel
character(len=40) :: Titorb
character(len=128) :: OrbFile
character(len=4000) :: BsLbl
character(len=3) :: PLab
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), allocatable :: DTmp(:)

!-- Banner.

ireturn = 99

!-- Define defaults and initialize.

call Init_ave(Title,iPrint,Wset,Wsum,PrOcc,PrEne,DensityBased,ThrOcc,Dummy(1),iDummy(1,1))

!-- Read input.

call Get_Averd_input(Title,Wset,iPrint,Nset,DensityBased,ThrOcc)

!-- Read some information from RUNFILE.

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
itBas = 0
do iSym=1,nSym
  itBas = itBas+nBas(isym)
end do
call Get_cArray('Unique Basis Names',BsLbl,(LENIN8)*itBas)

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

call GetMem('Overlap','Allo','Real',iS,lsmat+4)
OLabel = 'Mltpl  0'
irc = 0
iopt = 6
icomp = 1
isyml = 1
call RdOne(irc,iopt,OLabel,icomp,Work(iS),isyml)
if (iprint >= 99) then
  ind = 0
  do iSym=1,nSym
    call TriPrt('Overlap Matrix',' ',Work(iS+ind),nBas(iSym))
    ind = ind+nBas(iSym)*(nBas(iSym)+1)/2
  end do
end if

!-- Normalize weights.

do iset=1,mxsets
  Wsum = Wsum+wset(iset)
end do
do iset=1,Nset
  Wset(iset) = Wset(iset)/Wsum
end do

!-- Print some Bla Bla...

if (iprint >= 2) then
  call Print_Input(Title,nSym,nBas,wSet,nSet)
end if

!-- Do the dirty work. Different paths for orbital- and density-based
!   averageing.

call GetMem('Density','Allo','Real',iDao,ntot2)
call dcopy_(ntot2,[Zero],0,Work(iDao),1)
if (.not. DensityBased) then
  Luinp = 10
  call GetMem('Orbitals','Allo','Real',iCMO,ntot2)
  call GetMem('Occ','Allo','Real',iOcc,ntot)
  do iset=1,Nset
    Fname = 'NAT001'
    write(Fname(4:6),'(i3.3)') iset
    ! Read orbital coefficients and occupation numbers.
    call RdVec(Fname,Luinp,'CO',Nsym,nBas,nBas,Work(iCMO),Work(iOcc),Dummy,iDummy,Titorb,0,iErr)
    iC = 0
    iO = 0
    iD = 0
    ! Up-date average density matrix.
    do isym=1,nSym
      kaunter = 0
      do i=1,nBas(iSym)
        do j=1,nBas(iSym)
          do k=1,nBas(iSym)
            Work(iDao+iD+kaunter) = Work(iDao+iD+kaunter)+Wset(iSet)*Work(iOcc+iO+k-1)*Work(iCMO+iC+i+(k-1)*nBas(iSym)-1)* &
                                    Work(iCMO+iC+j+(k-1)*nBas(iSym)-1)
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
      ThrO = 1d-5
      call Primo(Titorb,PrOcc,PrEne,ThrO,Dummy(1),nSym,nBas,nBas,BsLbl,Dummy,Work(iOcc),Work(iCMO),-1)
    end if
  end do
  call GetMem('Orbitals','Free','Real',iCMO,ntot2)
  call GetMem('Occ','Free','Real',iOcc,ntot)
else
  call GetMem('DensityT','Allo','Real',iDtemp,lsmat)
  call dcopy_(lsmat,[Zero],0,Work(iDtemp),1)
  do iset=1,Nset
    Fname = 'RUN001'
    write(Fname(4:6),'(i3.3)') iset
    call NameRun(Fname)
    ! Collect density from runfile.
    call mma_allocate(DTmp,lsmat,Label='DTmp')
    call Get_D1ao(Dtmp,lsmat)
    call DaxPy_(lsmat,Wset(iset),Dtmp,1,Work(iDtemp),1)
    call mma_deallocate(DTmp)
  end do
  ! Square the density matrix.
  iDt = 0
  iDs = 0
  do iSym=1,nSym
    nB = nBas(iSym)
    call Dsq(Work(iDtemp+iDt),Work(iDao+iDs),1,nB,nB)
    iDt = iDt+nB*(nB+1)/2
    iDs = iDs+nB**2
  end do
  call GetMem('DensityT','Free','Real',iDtemp,lsmat)
end if

!-- With the average density in store, lets orthogonalize (canonical),
!   then diagonalize to get natural orbitals.

indT = 0
indS = 0
indB = 0
call GetMem('NatOrbAcc','Allo','Real',iOrbs,ntot2)
call GetMem('NatOccAcc','Allo','Real',iOccs,ntot)
do iSym=1,nSym
  nBT = nBas(iSym)*(nBas(iSym)+1)/2
  nBS = nBas(iSym)**2
  call GetMem('EigV','Allo','Real',iVecs,nBS)
  call GetMem('St','Allo','Real',iSt,nBT)
  call GetMem('Si','Allo','Real',iSi,nBT)
  call GetMem('Ss','Allo','Real',iSs,nBS)
  call GetMem('Sp','Allo','Real',iSp,nBS)
  call GetMem('AUX','Allo','Real',iAUX,nBS)
  call GetMem('TransS','Allo','Real',iTrans,nBS)
  call GetMem('TransSi','Allo','Real',iTrani,nBS)
  call GetMem('OrthoDensS','Allo','Real',iOrtoD,nBS)
  call GetMem('OrthoDensT','Allo','Real',iOrtoDt,nBT)
  call GetMem('Occs','Allo','Real',iOccNat,nBas(iSym))
  call dcopy_(nBT,[Zero],0,Work(iSt),1)
  call dcopy_(nBT,[Zero],0,Work(iSi),1)
  kaunter = 0
  do iB1=1,nBas(iSym)
    do iB2=1,nBas(iSym)
      Work(iVecs+kaunter) = Zero
      if (iB1 == iB2) Work(iVecs+kaunter) = One
      kaunter = kaunter+1
    end do
  end do
  call Jacob(Work(iS+indT),Work(iVecs),nBas(iSym),nBas(iSym))
  do i=1,nBas(iSym)
    Sqroot = sqrt(Work(iS+indT+i*(i+1)/2-1))
    Work(iSt+i*(i+1)/2-1) = Sqroot
    Work(iSi+i*(i+1)/2-1) = One/Sqroot
  end do
  call Square(Work(iSt),Work(iSs),1,nBas(iSym),nBas(iSym))
  call Square(Work(iSi),Work(iSp),1,nBas(iSym),nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iVecs),nBas(iSym),Work(iSs),nBas(iSym),Zero,Work(iAUX),nBas(iSym))
  call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iAUX),nBas(iSym),Work(iVecs),nBas(iSym),Zero,Work(iTrans), &
              nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iVecs),nBas(iSym),Work(iSp),nBas(iSym),Zero,Work(iAUX),nBas(iSym))
  call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iAUX),nBas(iSym),Work(iVecs),nBas(iSym),Zero,Work(iTrani), &
              nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iTrans),nBas(iSym),Work(iDao+indS),nBas(iSym),Zero,Work(iAUX), &
              nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iAUX),nBas(iSym),Work(iTrans),nBas(iSym),Zero,Work(iOrtoD), &
              nBas(iSym))
  kaunter = 0
  do iB1=1,nBas(iSym)
    do iB2=1,nBas(iSym)
      Work(iVecs+kaunter) = Zero
      if (iB1 == iB2) Work(iVecs+kaunter) = One
      kaunter = kaunter+1
    end do
  end do
  kaunter = 0
  do i=1,nBas(iSym)
    do j=1,i
      Work(iOrtoDt+kaunter) = Work(iOrtoD+i+(j-1)*nBas(iSym)-1)
      kaunter = kaunter+1
    end do
  end do
  call Jacob(Work(iOrtoDt),Work(iVecs),nBas(iSym),nBas(iSym))
  call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Work(iTrani),nBas(iSym),Work(iVecs),nBas(iSym),Zero,Work(iAUX), &
              nBas(iSym))
  kaunt = 0
  kaunter = 0
  do i=1,nBas(iSym)
    do j=1,i
      if (i == j) then
        Work(iOccNat+kaunt) = Work(iOrtoDt+kaunter)
        kaunt = kaunt+1
      end if
      kaunter = kaunter+1
    end do
  end do
  call Jacord3(Work(iOccNat),Work(iAUX),nBas(iSym),nBas(iSym))
  call Add_Info('AVERD_OCC',Work(iOccNat),5,5)
  if (iPrint >= 5) then
    write(Titorb,'(a)') 'All average orbitals in this irrep.'
    Thr = -1d0
    call Primo(Titorb,.true.,.false.,Thr,Dum,1,nBas(iSym),nBas(iSym),BsLbl,Dummy,Work(iOccNat),Work(iAUX),-1)
  end if
  call dcopy_(nBS,Work(iAUX),1,Work(iOrbs+indS),1)
  call dcopy_(nBas(iSym),Work(iOccNat),1,Work(iOccs+indB),1)
  call GetMem('EigV','Free','Real',iVecs,nBS)
  call GetMem('St','Free','Real',iSt,nBT)
  call GetMem('Si','Free','Real',iSi,nBT)
  call GetMem('Ss','Free','Real',iSs,nBS)
  call GetMem('Sp','Free','Real',iSp,nBS)
  call GetMem('AUX','Free','Real',iAUX,nBS)
  call GetMem('TransS','Free','Real',iTrans,nBS)
  call GetMem('TransSi','Free','Real',iTrani,nBS)
  call GetMem('OrthoDensS','Free','Real',iOrtoD,nBS)
  call GetMem('OrthoDensT','Free','Real',iOrtoDt,nBT)
  call GetMem('Occs','Free','Real',iOccNat,nBas(iSym))
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
call GetMem('Zeros','Allo','Real',iZero,ntot)
call dcopy_(ntot,[Zero],0,Work(iZero),1)
LuOut = 65
LuOut = IsFreeUnit(LuOut)
Title = 'Average Orbitals'
OrbFile = 'AVEORB'
Plab = 'COE'
call WrVec(OrbFile,LuOut,Plab,nSym,nBas,nBas,Work(iOrbs),Work(iOccs),Work(iZero),iDummy,Title)

!-- Say something about orbital occupation.

write(u6,*)
write(u6,*)
write(u6,'(a)') ' |  Average orbital occupation.'
write(u6,'(a)') ' |-----------------------------'
write(u6,'(a,e18.8)') ' |    Threshold: ',ThrOcc
write(u6,*)
nOrb = 0
iO = 0
do iSym=1,nSym
  do iB=1,nBas(iSym)
    if (Work(iOccs+iO+iB-1) < ThrOcc) cycle
    nOrb = nOrb+1
  end do
  write(u6,'(a,i2,a,i4)') '      Symmetry:',iSym,'   Number of orbitals below threshold:',nOrb
  iO = iO+nBas(iSym)
end do
write(u6,*)
call GetMem('Zeros','Free','Real',iZero,ntot)
call GetMem('NatOrbAcc','Free','Real',iOrbs,ntot2)
call GetMem('NatOccAcc','Free','Real',iOccs,ntot)
call GetMem('Density','Free','Real',iDao,ntot2)
call GetMem('Overlap','Free','Real',iS,lsmat+4)

!-- Good Bye.

ireturn = 0

return

end subroutine Averd
