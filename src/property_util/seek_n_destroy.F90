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

subroutine Seek_n_Destroy(nBasAtoms,nDimSubD,ipSubVal,ipSubVec,nBast,Threshold,ThrsD,TotElec,ipSubDNAOindex,ipDS,iWhat,ipAtomA, &
                          ipAtomB,iCounter,IAtom,JAtom,ipWhat,ipAtomC,KAtom)
!----------------------------------------------------------------------*
! Seek n Destroy: search the sub matrix for eigen values >= Thrs and   *
! deplete the corresponding scaled eigenvector from the main matrix.   *
!----------------------------------------------------------------------*

use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nBasAtoms, nDimSubD, ipSubVal, ipSubVec, nBast, ipSubDNAOindex, ipDS, iWhat, ipAtomA, ipAtomB, &
                                 IAtom, JAtom, ipWhat, ipAtomC, KAtom
real(kind=wp), intent(in) :: Threshold, ThrsD
real(kind=wp), intent(inout) :: TotElec
integer(kind=iwp), intent(inout) :: iCounter
integer(kind=iwp) :: I, iCounterOld, iCounterTrue, iFound, iFoundOrb, ipEiVal, ipGood, ipScrM, ipScrV, iStHas2bFnd, J, K
real(kind=wp) :: Accumulate, E, EigenNorm, Thrs, Thrs_Original, TotElecAvail, TotElecCount, TotElecFound
real(kind=wp), external :: DNRM2_
#include "WrkSpc.fh"

Thrs = Threshold
Thrs_Original = Threshold

!----------------------------------------------------------------------*
! Now loop over the eigenvalues looking for orbitals >= Thrs

call Allocate_iWork(ipGood,nBasAtoms)
call Allocate_Work(ipEiVal,nBasAtoms)
call IZero(iWork(ipGood),nBasAtoms)
call FZero(Work(ipEiVal),nBasAtoms)

iFoundOrb = 0
iStHas2bFnd = 0
TotElecCount = Zero
TotElecAvail = Zero
TotElecFound = Zero

do I=1,nBasAtoms
  TotElecAvail = TotElecAvail+Work(ipSubVal+I-1)
end do

if ((TotElecAvail > 2*ThrsD) .and. (iWhat == 1) .and. (Thrs >= 1.999_wp)) iStHas2bFnd = 1

if (TotElecAvail < Zero) goto 666

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Total electrons in eigenvalues = ',TotElecAvail
write(u6,*) 'Something has to be found = ',iStHas2bFnd
#endif

E = Zero
do K=1,NBAST
  E = E+Work(ipDS+(K-1)*NBAST+K-1)
end do

if (E < Zero) goto 666

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#endif

!----------------------------------------------------------------------*
777 continue
!----------------------------------------------------------------------*

do I=1,nBasAtoms
  if ((Work(ipSubVal+I-1) > Thrs) .and. (Work(ipSubVal+I-1) <= ThrsD)) then

    ! One good orbital found
    ! iFoundOrb: number of Good orbitals found
    ! iWork(ipGood): which of the eigenvectors corresponds to the big eigenvalue
    ! Work(ipEiVal): the big eigenvalue

    iFoundOrb = iFoundOrb+1
    iWork(ipGood+iFoundOrb-1) = I
    Work(ipEiVal+iFoundOrb-1) = Work(ipSubVal+I-1)

    ! It shouldn't be necessary, but in some cases (like radicals) a bit
    ! more than 2 electrons are found in core or lone pair orbitals. This
    ! way everything comes out cleaner. Possibly I'll try to remove this
    ! later.

    if (Work(ipEiVal+iFoundOrb-1) > 2) Work(ipEiVal+iFoundOrb-1) = 2

    TotElecFound = TotElecFound+Work(ipEiVal+iFoundOrb-1)
  end if
end do

#ifdef _DEBUGPRINT_
if (TotElecFound > TotElecAvail) then
  write(u6,*)
  write(u6,*) 'Something fishy is going on'
  write(u6,*) 'iFoundOrb     = ',iFoundOrb
  write(u6,*) 'TotElecFound  = ',TotElecFound
  write(u6,*) 'TotElecAvail  = ',TotElecAvail
end if
#endif

! If we didn't find any bond, and something had to be found anyway,
! we give it another try, with lower threshold, just once

if ((iFoundOrb == 0) .and. (iStHas2bFnd == 1)) then
  iStHas2bFnd = 0
  !Thrs = Thrs-0.05_wp
  !Thrs = Thrs-0.01_wp
  Thrs = Thrs-0.005_wp
  goto 777
end if

Thrs = Thrs_Original

TotElecCount = TotElecFound-TotElecAvail

if ((iFoundOrb > 0) .and. (TotElecCount <= 0.1_wp)) then

  TotElec = TotElec+TotElecFound

  if (iWhat == 1) then
    iCounterOld = iCounter
    iCounterTrue = iCounter

    if (iCounter > 0) then
      iFound = -1
      do I=0,iCounter-1
        if ((iWork(ipAtomA+I) == IAtom) .and. (iWork(ipAtomB+I) == JAtom)) then
          iFound = I
        end if
      end do
      if (iFound >= 0) then
        iCounterTrue = iFound
        iCounterOld = iCounterOld-1
      end if
    end if

    iWork(ipAtomA+iCounterTrue) = IAtom
    iWork(ipAtomB+iCounterTrue) = JAtom
    Accumulate = Zero
    do I=1,iFoundOrb
      Accumulate = Accumulate+Work(ipEiVal+I-1)
    end do
    Work(ipWhat+iCounterTrue) = Work(ipWhat+iCounterTrue)+(Accumulate/2)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'NBO bond order between Atom ',iWork(ipAtomA+iCounterTrue),' and Atom ',iWork(ipAtomB+iCounterTrue),' is ', &
                Work(ipWhat+iCounterTrue)
#   endif

    iCounter = iCounterOld+1
  end if

  if (iWhat == 2) then
    do I=1,iFoundOrb
      iWork(ipAtomA+iCounter) = IAtom
      Work(ipWhat+iCounter) = Work(ipEiVal+I-1)

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'NBO located ',Work(ipEiVal+I-1),' electrons'
      write(u6,*) 'non bonding on atom ',iWork(ipAtomA+iCounter)
#     endif

      iCounter = iCounter+1
    end do
  end if

  if (iWhat == 3) then
    do I=1,iFoundOrb
      iWork(ipAtomA+iCounter) = IAtom
      iWork(ipAtomB+iCounter) = JAtom
      iWork(ipAtomC+iCounter) = KAtom
      Work(ipWhat+iCounter) = Work(ipEiVal+I-1)/2

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'NBO tricenter bond order between Atom ',iWork(ipAtomA+iCounter),' and Atom ',iWork(ipAtomB+iCounter), &
                  ' and Atom ',iWork(ipAtomC+iCounter),' is ',Work(ipWhat+iCounter)
#     endif

      iCounter = iCounter+1
    end do
  end if

  do I=1,iFoundOrb
    call Allocate_Work(ipScrM,nDimSubD)
    call Allocate_Work(ipScrV,nBasAtoms)

    call FZero(Work(ipScrM),nDimSubD)
    call FZero(Work(ipScrV),nBasAtoms)

    do J=1,nBasAtoms
      Work(ipScrV+J-1) = Work(ipSubVec+(iWork(ipGood+I-1)-1)*nBasAtoms+J-1)
    end do

#   ifdef _DEBUGPRINT_
    call RecPrt('Good orbital eigenvector',' ',Work(ipScrV),nBasAtoms,1)
#   endif

    EigenNorm = DNRM2_(nBasAtoms,Work(ipScrV),1)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Vector euclidean norm = ',EigenNorm
#   endif

    call DScal_(nBasAtoms,1/EigenNorm,Work(ipScrV),1)

#   ifdef _DEBUGPRINT_
    write(u6,*)
    call RecPrt('Normalized orbital eigenvector',' ',Work(ipScrV),nBasAtoms,1)
#   endif

    call DGEMM_('N','T',nBasAtoms,nBasAtoms,1,One,Work(ipScrV),nBasAtoms,Work(ipScrV),nBasAtoms,Zero,Work(ipScrM),nBasAtoms)

#   ifdef _DEBUGPRINT_
    call RecPrt('Multiplied Eigenvectors',' ',Work(ipScrM),nBasAtoms,nBasAtoms)

    write(u6,*) 'Eigen value to be multiplied = ',Work(ipEiVal+I-1)
#   endif

    call DScal_(nDimSubD,Work(ipEiVal+I-1),Work(ipScrM),1)

#   ifdef _DEBUGPRINT_
    call RecPrt('Scaled Mult eigenvector',' ',Work(ipScrM),nBasAtoms,nBasAtoms)
#   endif

    ! DS Matrix depletion not just along the diagonal

    do J=0,nDimSubD-1
      K = iWork(ipSubDNAOindex+J)
      Work(ipDS+K) = Work(ipDS+K)-Work(ipScrM+J)
    end do

#   ifdef _DEBUGPRINT_
    call RecPrt('DS-NAO depleted Matrix = ',' ',Work(ipDS),NBAST,NBAST)
    E = Zero
    do K=1,NBAST
      E = E+Work(ipDS+(K-1)*NBAST+K-1)
    end do
    write(u6,*)
    write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#   endif

    call Free_Work(ipScrV)
    call Free_Work(ipScrM)
  end do
end if

666 continue
call Free_iWork(ipGood)
call Free_Work(ipEiVal)

return

end subroutine Seek_n_Destroy
