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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine DM_FNO_RHF(irc,nSym,nBas,nFro,nIsh,nSsh,nDel,CMOI,EOcc,EVir,DM0,DM)
!***********************************************************************
!                                                                      *
!     Purpose:  setup of FNO density matrix calculation (RHF-based)    *
!                                                                      *
!     Author:   F. Aquilante  (Geneva, Sep 2010)                       *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use ChoMP2, only: MP2_small
use Molcas, only: MxBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: iRC
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nSsh(nSym), nDel(nSym)
real(kind=wp), intent(in) :: CMOI(*), EOcc(*), EVir(*)
real(kind=wp), intent(_OUT_) :: DM0(*), DM(*)
integer(kind=iwp) :: i, ifr, ioff, iSkip, iSym, iTo, j, jD, jOcc, jOff, jp, jTo, jVir, kDM, kfr, kij, kOff, kTo, lij, lnDel(8), &
                     lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), lOff, nBasT, nBmx, nCMO, nOA, nOkk, nOrb, nSQ, nTri, nVV
real(kind=wp) :: Dummy, SqOcc, tmp
real(kind=wp), allocatable :: CMO(:,:), DMAT(:), EOrb(:,:)

irc = 0
MP2_small = .false.

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nBasT = 0
ntri = 0
nSQ = 0
nBmx = 0
nOrb = 0
nVV = 0
do i=1,nSym
  nBasT = nBasT+nBas(i)
  nOrb = nOrb+nFro(i)+nIsh(i)+nSsh(i)+nDel(i)
  ntri = ntri+nTri_Elem(nBas(i))
  nSQ = nSQ+nBas(i)**2
  nVV = nVV+nSsh(i)**2
  nBmx = max(nBmx,nBas(i))
end do
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

NCMO = nSQ
call mma_allocate(CMO,nCMO,2,Label='CMO')
CMO(:,1) = CMOI(1:nCMO)

nOA = 0
do iSym=1,nSym  ! setup info
  lnFro(iSym) = nFro(iSym)
  lnOcc(iSym) = nIsh(iSym)
  nOA = nOA+lnOcc(iSym)
  lnVir(iSym) = nSsh(iSym)
  lnOrb(iSym) = lnOcc(iSym)+lnVir(iSym)
  lnDel(iSym) = nDel(iSym)
end do

call mma_Allocate(EOrb,nOrb,4,Label='EOrb')
jOff = 0
kOff = 0
lOff = 0
do iSym=1,nSym
  jp = 1+lOff+nFro(iSym)
  jOcc = jOff+1
  EOrb(jp:jp+nIsh(iSym)-1,1) = EOcc(jOcc:jOcc+nIsh(iSym)-1)
  jVir = kOff+1
  jp = jp+nIsh(iSym)
  EOrb(jp:jp+nSsh(iSym)-1,1) = EVir(jVir:jVir+nSsh(iSym)-1)
  jOff = jOff+nIsh(iSym)
  kOff = kOff+nSsh(iSym)
  lOff = lOff+nBas(iSym)
end do
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = 1+ioff+nFro(iSym)
  ito = 1+joff
  EOrb(ito:ito+nIsh(iSym)-1,3) = EOrb(ifr:ifr+nIsh(iSym)-1,1)
  ifr = 1+ioff+nFro(iSym)+nIsh(iSym)
  ito = 1+koff
  EOrb(ito:ito+nSsh(iSym)-1,4) = EOrb(ifr:ifr+nSsh(iSym)-1,1)
  ioff = ioff+nBas(iSym)
  joff = joff+nIsh(iSym)
  koff = koff+nSsh(iSym)
end do
call mma_Allocate(DMAT,nVV+nOA,Label='DMAT')
DMAT(:) = Zero

call FnoSCF_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)

CMO(:,2) = Zero
iOff = 0
do iSym=1,nSym
  kfr = 1+iOff+nBas(iSym)*nFro(iSym)
  kto = 1+iOff+nBas(iSym)*lnFro(iSym)
  CMO(kto:kto+nBas(iSym)*lnOcc(iSym)-1,2) = CMO(kfr:kfr+nBas(iSym)*lnOcc(iSym)-1,1)
  kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  CMO(kto:kto+nBas(iSym)*lnVir(iSym)-1,2) = CMO(kfr:kfr+nBas(iSym)*lnVir(iSym)-1,1)
  iOff = iOff+nBas(iSym)**2
end do
call Check_Amp_SCF(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,CMO(:,2),EOrb(:,3),EOrb(:,4),DMAT(1:nVV),DMAT(nVV+1:))
  if (irc /= 0) then
    write(u6,*) 'MP2 pseudodensity calculation failed !'
    call Abend()
  end if
else
  write(u6,*)
  write(u6,*) 'There are ZERO amplitudes T(ai,bj) with the given '
  write(u6,*) 'combinations of inactive and virtual orbitals !! '
  write(u6,*) 'Check your input and rerun the calculation! Bye!!'
  call Abend()
end if

! Compute the correlated density in AO basis
! -------------------------------------------------------------
jOcc = 1+nVV
!write(u6,*) ' Occ    : ',(DMAT(jOcc+j),j=0,nOA-1)
!write(u6,*) ' Sum    : ',ddot_(nOA,One,0,DMAT(jOcc),1)
DMAT(jOcc:jOcc+nOA-1) = Two*DMAT(jOcc:jOcc+nOA-1)+One

iOff = 0
jOff = 0
kDM = 1
do iSym=1,nSym

  kto = 1+jOff
  nOkk = nFro(iSym)+nIsh(iSym)
  call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk, &
                 Two,CMO(kto,1),nBas(iSym), &
                 CMO(kto,1),nBas(iSym), &
                 Zero,DM0(kDM),nBas(iSym))

  sqocc = sqrt(Two)
  CMO(kto:kto+nBas(iSym)*nFro(iSym)-1,1) = sqocc*CMO(kto:kto+nBas(iSym)*nFro(iSym)-1,1)
  do j=0,nIsh(iSym)-1
    sqocc = sqrt(DMAT(jOcc+j))
    ito = kto+nBas(iSym)*j
    CMO(ito:ito+nBas(iSym)-1,1) = sqocc*CMO(ito:ito+nBas(iSym)-1,1)
  end do
  call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk, &
                 One,CMO(kto,1),nBas(iSym), &
                 CMO(kto,1),nBas(iSym), &
                 Zero,DM(kDM),nBas(iSym))

  if (nSsh(iSym) > 0) then
    jD = 1+iOff
    ! Eigenvectors will be in increasing order of eigenvalues
    call Eigen_Molcas(nSsh(iSym),DMAT(jD),EOrb(:,2),Eorb(:,1))
    ! Reorder to get relevant eigenpairs first
    do j=1,nSsh(iSym)/2
      do i=1,nSsh(iSym)
        lij = jD-1+nSsh(iSym)*(j-1)+i
        kij = jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
        tmp = DMAT(lij)
        DMAT(lij) = DMAT(kij)
        DMAT(kij) = tmp
      end do
      tmp = EOrb(j,2)
      EOrb(j,2) = EOrb(nSsh(iSym)-j,2)
      EOrb(nSsh(iSym)-j,2) = tmp
    end do

    ! Compute new MO coeff. : X=C*U
    kfr = 1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    kto = 1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym), &
                One,CMO(kfr,2),nBas(iSym), &
                DMAT(jD),nSsh(iSym), &
                Zero,CMO(kto,1),nBas(iSym))

    !write(u6,*) ' Occ_vir: ',(EOrb(j,2),j=1,nSsh(iSym))
    !write(u6,*) ' Sum_vir: ',ddot_(nSsh(iSym),One,0,EOrb(:,2),1)
    do j=0,nSsh(iSym)-1
      sqocc = sqrt(Two*EOrb(1+j,2))
      jto = kto+nBas(iSym)*j
      CMO(jto:jto+nBas(iSym)-1,1) = sqocc*CMO(jto:jto+nBas(iSym)-1,1)
    end do
    call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nSsh(iSym), &
                   One,CMO(kto,1),nBas(iSym), &
                   CMO(kto,1),nBas(iSym), &
                   One,DM(kDM),nBas(iSym))

    iOff = iOff+nSsh(iSym)**2
  end if
  jOff = jOff+nBas(iSym)**2
  kDM = kDM+nTri_Elem(nBas(iSym))
  jOcc = jOcc+nIsh(iSym)
end do

call mma_deAllocate(EOrb)
call mma_deAllocate(DMAT)
call mma_deallocate(CMO)

return

end subroutine DM_FNO_RHF
