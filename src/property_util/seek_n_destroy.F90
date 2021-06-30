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
! Copyright (C) 2006, Luca De Vico                                     *
!               Valera Veryazov                                        *
!***********************************************************************

subroutine Seek_n_Destroy(nBasAtoms,SubVal,SubVec,nBast,Threshold,ThrsD,TotElec,SubDNAOindex,DS,iWhat,AtomA,AtomB,iCounter,IAtom, &
                          JAtom,What,AtomC,KAtom)
!----------------------------------------------------------------------*
! Seek n Destroy: search the sub matrix for eigen values >= Thrs and   *
! deplete the corresponding scaled eigenvector from the main matrix.   *
!----------------------------------------------------------------------*

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nBasAtoms, nBast, SubDNAOindex(2,nBasAtoms,nBasAtoms), iWhat, IAtom, JAtom, KAtom
real(kind=wp), intent(in) :: SubVal(nBasAtoms), SubVec(nBasAtoms,nBasAtoms), Threshold, ThrsD
real(kind=wp), intent(inout) :: TotElec, DS(nBasAtoms,nBasAtoms), What(*)
integer(kind=iwp), intent(inout) :: AtomA(*), AtomB(*), iCounter, AtomC(*)
integer(kind=iwp) :: I, iCounterOld, iCounterTrue, iFound, iFoundOrb, iStHas2bFnd, J, JJ, K, KK
real(kind=wp) :: Accumulate, E, EigenNorm, Thrs, Thrs_Original, TotElecAvail, TotElecCount, TotElecFound
integer(kind=iwp), allocatable :: Good(:)
real(kind=wp), allocatable :: EiVal(:), ScrM(:,:), ScrV(:)
real(kind=wp), external :: DNRM2_

Thrs = Threshold
Thrs_Original = Threshold

!----------------------------------------------------------------------*
! Now loop over the eigenvalues looking for orbitals >= Thrs

call mma_allocate(Good,nBasAtoms,label='Good')
call mma_allocate(EiVal,nBasAtoms,label='EiVal')
Good(:) = 0
EiVal(:) = Zero

iFoundOrb = 0
iStHas2bFnd = 0
TotElecCount = Zero
TotElecAvail = Zero
TotElecFound = Zero

do I=1,nBasAtoms
  TotElecAvail = TotElecAvail+SubVal(I)
end do

if ((TotElecAvail > 2*ThrsD) .and. (iWhat == 1) .and. (Thrs >= 1.999_wp)) iStHas2bFnd = 1

if (TotElecAvail < Zero) then
  call End1()
  return
end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Total electrons in eigenvalues = ',TotElecAvail
write(u6,*) 'Something has to be found = ',iStHas2bFnd
#endif

E = Zero
do K=1,NBAST
  E = E+DS(K,K)
end do

if (E < Zero) then
  call End1()
  return
end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#endif

!----------------------------------------------------------------------*
do

  do I=1,nBasAtoms
    if ((SubVal(I) > Thrs) .and. (SubVal(I) <= ThrsD)) then

      ! One good orbital found
      ! iFoundOrb: number of Good orbitals found
      ! Good: which of the eigenvectors corresponds to the big eigenvalue
      ! EiVal: the big eigenvalue

      iFoundOrb = iFoundOrb+1
      Good(iFoundOrb) = I
      EiVal(iFoundOrb) = SubVal(I)

      ! It shouldn't be necessary, but in some cases (like radicals) a bit
      ! more than 2 electrons are found in core or lone pair orbitals. This
      ! way everything comes out cleaner. Possibly I'll try to remove this
      ! later.

      if (EiVal(iFoundOrb) > 2) EiVal(iFoundOrb) = 2

      TotElecFound = TotElecFound+EiVal(iFoundOrb)
    end if
  end do

# ifdef _DEBUGPRINT_
  if (TotElecFound > TotElecAvail) then
    write(u6,*)
    write(u6,*) 'Something fishy is going on'
    write(u6,*) 'iFoundOrb     = ',iFoundOrb
    write(u6,*) 'TotElecFound  = ',TotElecFound
    write(u6,*) 'TotElecAvail  = ',TotElecAvail
  end if
# endif

  ! If we didn't find any bond, and something had to be found anyway,
  ! we give it another try, with lower threshold, just once

  if ((iFoundOrb == 0) .and. (iStHas2bFnd == 1)) then
    iStHas2bFnd = 0
    !Thrs = Thrs-0.05_wp
    !Thrs = Thrs-0.01_wp
    Thrs = Thrs-0.005_wp
  else
    exit
  end if
end do

Thrs = Thrs_Original

TotElecCount = TotElecFound-TotElecAvail

if ((iFoundOrb > 0) .and. (TotElecCount <= 0.1_wp)) then

  TotElec = TotElec+TotElecFound

  if (iWhat == 1) then
    iCounterOld = iCounter
    iCounterTrue = iCounter+1

    if (iCounter > 0) then
      iFound = -1
      do I=1,iCounter
        if ((AtomA(I) == IAtom) .and. (AtomB(I) == JAtom)) then
          iFound = I
        end if
      end do
      if (iFound >= 0) then
        iCounterTrue = iFound
        iCounterOld = iCounterOld-1
      end if
    end if

    AtomA(iCounterTrue) = IAtom
    AtomB(iCounterTrue) = JAtom
    Accumulate = Zero
    do I=1,iFoundOrb
      Accumulate = Accumulate+EiVal(I)
    end do
    What(iCounterTrue) = What(iCounterTrue)+(Accumulate/2)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'NBO bond order between Atom ',AtomA(iCounterTrue),' and Atom ',AtomB(iCounterTrue),' is ',What(iCounterTrue)
#   endif

    iCounter = iCounterOld+1
  end if

  if (iWhat == 2) then
    do I=1,iFoundOrb
      iCounter = iCounter+1

      AtomA(iCounter) = IAtom
      What(iCounter) = EiVal(I)

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'NBO located ',EiVal(I),' electrons'
      write(u6,*) 'non bonding on atom ',AtomA(iCounter)
#     endif
    end do
  end if

  if (iWhat == 3) then
    do I=1,iFoundOrb
      iCounter = iCounter+1

      AtomA(iCounter) = IAtom
      AtomB(iCounter) = JAtom
      AtomC(iCounter) = KAtom
      What(iCounter) = EiVal(I)/2

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'NBO tricenter bond order between Atom ',AtomA(iCounter),' and Atom ',AtomB(iCounter),' and Atom ', &
                  AtomC(iCounter),' is ',What(iCounter)
#     endif
    end do
  end if

  do I=1,iFoundOrb
    call mma_allocate(ScrM,nBasAtoms,nBasAtoms,label='ScrM')
    call mma_allocate(ScrV,nBasAtoms,label='ScrV')

    ScrM(:,:) = Zero

    ScrV(:) = SubVec(:,Good(I))

#   ifdef _DEBUGPRINT_
    call RecPrt('Good orbital eigenvector',' ',ScrV,nBasAtoms,1)
#   endif

    EigenNorm = DNRM2_(nBasAtoms,ScrV,1)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Vector euclidean norm = ',EigenNorm
#   endif

    call DScal_(nBasAtoms,1/EigenNorm,ScrV,1)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    call RecPrt('Normalized orbital eigenvector',' ',ScrV,nBasAtoms,1)
#   endif

    call DGEMM_('N','T',nBasAtoms,nBasAtoms,1,One,ScrV,nBasAtoms,ScrV,nBasAtoms,Zero,ScrM,nBasAtoms)

#   ifdef _DEBUGPRINT_
    call RecPrt('Multiplied Eigenvectors',' ',ScrM,nBasAtoms,nBasAtoms)

    write(u6,*) 'Eigen value to be multiplied = ',EiVal(I)
#   endif

    call DScal_(nBasAtoms*nBasAtoms,EiVal(I),ScrM,1)

#   ifdef _DEBUGPRINT_
    call RecPrt('Scaled Mult eigenvector',' ',ScrM,nBasAtoms,nBasAtoms)
#   endif

    ! DS Matrix depletion not just along the diagonal

    do J=1,nBasAtoms
      do K=1,nBasAtoms
        KK = SubDNAOindex(1,K,J)
        JJ = SubDNAOindex(2,K,J)
        DS(KK,JJ) = DS(KK,JJ)-ScrM(K,J)
      end do
    end do

#   ifdef _DEBUGPRINT_
    call RecPrt('DS-NAO depleted Matrix = ',' ',DS,NBAST,NBAST)
    E = Zero
    do K=1,NBAST
      E = E+DS(K,K)
    end do
    write(u6,*)
    write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#   endif

    call mma_deallocate(ScrV)
    call mma_deallocate(ScrM)
  end do
end if

call End1()

return

contains

subroutine End1()
  call mma_deallocate(Good)
  call mma_deallocate(EiVal)
end subroutine End1

end subroutine Seek_n_Destroy
