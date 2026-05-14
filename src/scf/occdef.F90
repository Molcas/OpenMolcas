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

!#define _DEBUGPRINT_
!#define _SPECIAL_DEBUGPRINT_
subroutine OccDef(Occ,mmB,nD,CMO,mBB)

#include "compiler_features.h"

use InfSCF, only: kOV, mOV, nBas, nFro, nnb, nOcc, nOrb, nOV, nSym, OccSet_e, OccSet_m, OnlyProp, OrbType
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: mmB, nD, mBB
real(kind=wp), intent(out) :: Occ(mmB,nD)
real(kind=wp), target, intent(inout) :: CMO(mBB,nD)
integer(kind=iwp) :: iD, iOcc, iOff, iOrb, iSym, iTmp, jEOr, jOff, jOrb, MaxnOcc, MinnOcc, mOcc, mSet, Muon_i, nB, nOcc_e, nOcc_m
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: Tmp
integer(kind=iwp), allocatable :: iFerm(:)
real(kind=wp), allocatable :: OccTmp(:)
real(kind=wp), pointer :: pCMO(:,:)

! Form occupation numbers
!
! Note, this puts only in the occupation numbers of the electronic
! orbitals. The muonic occupation numbers are put in place in
! NewOrb_Scf on the first iteration.

if (OnlyProp) return

#ifdef _DEBUGPRINT_
do iD=1,nD
  write(u6,*) 'iD=',iD
  write(u6,'(a,8i5)') 'sorb: nOcc   ',(nOcc(i,iD),i=1,nSym)
end do
if (allocated(OccSet_e)) then
  do iD=1,nD
    write(u6,*) 'iD=',iD
    iOff = 1
    do iSym=1,nSym
      call RecPrt('OccSet_e',' ',OccSet_e(iOff,iD),1,nOcc(iSym,iD))
      iOff = iOff+nOcc(iSym,iD)
    end do
  end do
end if
if (allocated(OccSet_m)) then
  do iD=1,nD
    write(u6,*) 'iD=',iD
    iOff = 1
    do iSym=1,nSym
      call RecPrt('OccSet_m',' ',OccSet_m(iOff,iD),1,nOcc(iSym,iD))
      iOff = iOff+nOcc(iSym,iD)
    end do
  end do
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! First we set up the electronic occupation numbers either as
! by default or as specied by the user. The numbers are stored
! in the array Occ.

Occ(:,:) = Zero

if (.not. allocated(OccSet_e)) then

  ! Default is that occupied orbitals are set to 2 and 1
  ! for RHF and UHF, respectively.

  do iD=1,nD
    mOcc = 0
    do iSym=1,nSym
      Occ(mOcc+1:mOcc+nOcc(iSym,iD),iD) = real(3-nD,kind=wp) ! Value 2 or 1
      mOcc = mOcc+nOrb(iSym)
    end do
  end do

else

  ! User define occupation numbers according to the OCCN key word.

  do iD=1,nD
    mOcc = 0
    mSet = 0
    do iSym=1,nSym
      Occ(mOcc+1:mOcc+nOcc(iSym,iD),iD) = OccSet_e(mSet+1:mSet+nOcc(iSym,iD),iD)
      mOcc = mOcc+nOrb(iSym)
      mSet = mSet+nOcc(iSym,iD)
    end do
  end do
  call mma_deallocate(OccSet_e)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Here we put in the occupations from the muons

! Here we will have to run through the orbitals and identify
! if the orbitals are muonic or electronic. As we proceed we
! will fill up the vector with the occupation number from
! the list of either muonic or electronic occupation numbers.

if (allocated(OccSet_m)) then

  call mma_Allocate(OccTmp,mmB,Label='OccTmp')
  call mma_allocate(iFerm,nnB,Label='iFerm')
  call Get_iArray('Fermion IDs',iFerm,nnB)
  !write(u6,*) 'iFerm=',iFerm

  do iD=1,nD

    ! Store the electronic occupation numbers in OccTmp

    OccTmp(:) = Occ(:,iD)
    !write(u6,*) 'OccTmp=',OccTmp
    !write(u6,*) 'Occset_m=',Occset_m
    Occ(:,iD) = Zero

    nOcc_e = 0   ! number of occupied electronic orbitals
    nOcc_m = 0   ! number of occupied muonic orbitals
    iOff = 1     ! Offset to the symmetry sections of CMO
    jEOr = 0     ! Offset to the array with occupation numbers
    do iSym=1,nSym

      ! map the relevant portion of CMO onto pCMO

      nB = nBas(iSym)
      pCMO(1:nB,1:nB) => CMO(iOff:iOff+nB**2-1,iD)

      do iOrb=1,nOrb(iSym)  ! Loop over the orbitals

        ! Compute identifier which is different from zero
        ! if the orbital is muonic.

        tmp = Zero
        tmp = sum(real(iFerm(jEOr+1:jEOr+nB),kind=wp)*abs(pCMO(1:nB,iOrb)))
        Muon_i = 0                  ! electronic
        if (tmp /= Zero) Muon_i = 1 ! muonic

        if (Muon_i == 0) then

          ! The orbital is electronic, put in an
          ! electronic occupation number. If the index
          ! is beyond the number of occupied orbitals
          ! put in a Zero.

          nOcc_e = nOcc_e+1
          if (nOcc_e <= nOcc(iSym,iD)) then
            Occ(jEor+iOrb,iD) = OccTmp(jEOr+nOcc_e)
          else
            Occ(jEor+iOrb,iD) = Zero
          end if
          !write(u6,*) 'Electronic:',iOrb,Occ(jEOr+iOrb,iD)

        else if (Muon_i == 1) then

          ! The orbital is muonic, put in an
          ! muonic occupation number. If the index
          ! is beyond the number of occupied orbitals
          ! put in a Zero.

          nOcc_m = nOcc_m+1
          if (nOcc_m <= nOcc(iSym,id)) then
            Occ(jEor+iOrb,iD) = OccSet_m(jEOr+nOcc_m,iD)
          else
            Occ(jEor+iOrb,iD) = Zero
          end if
          OrbType(jEor+iOrb,iD) = 1
          !write(u6,*) 'Muonic:',iOrb,Occ(jEOr+iOrb,iD)

        end if

      end do

      nullify(pCMO)
      iOff = iOff+nB**2
      jEOr = jEOr+nOrb(iSym)
    end do ! iSym

  end do   ! iD
# ifdef _SPECIAL_DEBUGPRINT_
  call DebugCMO(CMO,mBB,nD,Occ,mmB,nBas,nOrb,nSym,iFerm,'@ the last position')
# endif
  call mma_deallocate(iFerm)
  call mma_deAllocate(OccTmp)
  call mma_deallocate(OccSet_m)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Finally we will have to resort the orbitals with respect to their
! occupation numbers such that the orbitals which are formally in
! the occupied space but which are empty are reassigned to being
! virtual orbitals.
!
! This means that the the CMO and Occ arrays will be resorted, and
! that the nOcc array will be updated with respect to the actual
! number of orbitals with occupation different from zero, i.e.
! they are virtual.

do iD=1,nD
# ifdef _DEBUGPRINT_
  write(u6,*) 'iD=',iD
  write(u6,*) 'nOccs(original):'
  write(u6,*) (nOcc(iSym,iD),iSym=1,nSym)
  write(u6,*) 'nOV=',nOV
# endif
  iOff = 1
  jOff = 0

  do iSym=1,nSym

    if (nOrb(iSym) == 0) cycle
    iOcc = 0
    do iOrb=1,nOrb(iSym)-1
      jOrb = iOrb+1
      if ((Occ(iOrb+jOff,iD) == Zero) .and. (Occ(jOrb+jOff,iD) > Occ(iOrb+jOff,iD))) then

        iTmp = OrbType(iOrb+jOff,iD)
        OrbType(iOrb+jOff,iD) = OrbType(jOrb+jOff,iD)
        OrbType(jOrb+jOff,iD) = iTmp

        Tmp = Occ(iOrb+jOff,iD)
        Occ(iOrb+jOff,iD) = Occ(jOrb+jOff,iD)
        Occ(jOrb+jOff,iD) = Tmp
        call DSwap_(nBas(iSym),CMO(iOff+(iOrb-1)*nBas(iSym),iD),1,CMO(iOff+(jOrb-1)*nBas(iSym),iD),1)
      end if
      if (Occ(iOrb+jOff,iD) /= Zero) iOcc = iOcc+1
    end do
    if (Occ(nOrb(iSym)+jOff,iD) /= Zero) iOcc = iOcc+1
    nOcc(iSym,iD) = iOcc

    jOff = jOff+nOrb(iSym)
    iOff = iOff+nBas(iSym)*nOrb(iSym)
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) 'iD=',iD
  write(u6,*) 'nOccs(new):'
  write(u6,*) (nOcc(iSym,iD),iSym=1,nSym)
# endif
end do

! Sort such that the occupied orbitals are first in each irrep.

call Sorb2CMOs(CMO,mBB,nD,Occ,mmB,nBas,nOrb,nSym,OrbType)
!                                                                      *
!***********************************************************************
!                                                                      *
! Recompute sizes since the nOcc array might have changed.

nOV = 0
mOV = 0
kOV(:) = 0
do iSym=1,nSym
  if (nD == 1) then
    maxnOcc = nOcc(iSym,1)
    minnOcc = nOcc(iSym,1)
  else
    maxnOcc = max(nOcc(iSym,1),nOcc(iSym,2))
    minnOcc = min(nOcc(iSym,1),nOcc(iSym,2))
  end if
  kOV(:) = kOV(:)+(nOcc(iSym,:)-nFro(iSym))*(nOrb(iSym)-nOcc(iSym,:))
  nOV = nOV+(maxnOcc-nFro(iSym))*(nOrb(iSym)-minnOcc)
end do
mOV = kOV(1)+kOV(2)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
do iD=1,nD
  iOff = 1
  jOff = 1
  do iSym=1,nSym
    call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
    call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),nBas(iSym),nOrb(iSym))
    do i=0,nOrb(iSym)-1
      write(u6,*) 'i,OrbType=',i+1,OrbType(iOff+i,iD)
    end do
    iOff = iOff+nOrb(iSym)
    jOff = jOff+nOrb(iSym)*nBas(iSym)
  end do
end do
#endif

return

#ifdef _SPECIAL_DEBUGPRINT_
contains

subroutine DebugCMO(CMO,nCMO,nD,Occ,nnB,nBas,nOrb,nSym,iFerm,Label)

  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: nCMO, nD, nnB, nSym, nBas(nSym), nOrb(nSym), iFerm(nnB)
  real(kind=wp), intent(in) :: CMO(nCMO,nD), Occ(nnB,nD)
  character(len=*) :: Label
  integer(kind=iwp) :: iD, iOff, iOrb, iSym, jOff
  real(kind=wp) :: tmp

  write(u6,*) Label
  do iD=1,nD
    write(u6,*)
    if (iD == 1) then
      if (nD == 1) then
        write(u6,*) ' RHF CMOs'
      else
        write(u6,*) ' UHF alpha CMOs'
      end if
    else
      write(u6,*) ' UHF beta CMOs'
    end if
    write(u6,*)
    jOff = 0
    iOff = 1
    do iSym=1,nSym
      do iOrb=1,nOrb(iSym)
        tmp = sum(real(iFerm(jOff+1:jOff+nBas(iSym)-1),kind=wp)*abs(CMO(iOff+(iOrb-1)*nBas(iSym):iOff+iOrb*nBas(iSym),iD))
        write(u6,*)
        if (tmp /= Zero) then
          write(u6,*) 'Muonic Orbital:',iOrb
        else
          write(u6,*) 'Electronic  Orbital:',iOrb
        end if
        write(u6,*) 'Occupation number:',Occ(jOff+iOrb,iD)
        call RecPrt('CMO',' ',CMO(iOff+(iOrb-1)*nBas(iSym),iD),1,nBas(iSym))
      end do
      jOff = jOff+nOrb(iSym)
      iOff = iOff+nBas(iSym)*nOrb(iSym)
    end do
  end do

  return

end subroutine DebugCMO
#endif

#undef _DEBUGPRINT_
end subroutine OccDef
