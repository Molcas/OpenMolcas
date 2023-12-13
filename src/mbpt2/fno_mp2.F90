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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine FNO_MP2(irc,nSym,nBas,nFro,nIsh,nSsh,nDel,CMOI,EOcc,EVir,vfrac,DoMP2,EMP2)
! Purpose:  setup of Frozen Natural Orbitals MP2 (FNO-MP2).
!
! Author:   F. Aquilante  (Geneva, Nov  2008)

use ChoMP2, only: DeMP2, MP2_small, shf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym)
integer(kind=iwp), intent(inout) :: nSsh(nSym), nDel(nSym)
real(kind=wp), intent(inout) :: CMOI(*), EVir(*)
real(kind=wp), intent(in) :: EOcc(*), vfrac
logical(kind=iwp), intent(in) :: DoMP2
real(kind=wp), intent(out) :: EMP2
integer(kind=iwp) :: i, ifr, ii, ioff, iSkip, iSym, ito, j, jD, jOcc, jOff, jp, jVir, k, kfr, kij, kOff, kto, lij, lnDel(8), &
                     lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), lOff, nAuxO(8), nBasT, nBmx, nBx, NCMO, nOA, nOrb, ns_V(8), nSQ, &
                     nSsh_t, nSx, ntri, nVV
real(kind=wp) :: Dummy, STrDF, STrDP, tmp, TrDF(8), TrDP(8)
integer(kind=iwp), allocatable :: iD(:)
real(kind=wp), allocatable :: CMO(:,:), OrbE(:,:), X(:)
real(kind=wp), external :: ddot_
#include "Molcas.fh"

irc = 0
MP2_small = .false.
shf = Zero

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
  ns_V(i) = 0
  nBasT = nBasT+nBas(i)
  nOrb = nOrb+nFro(i)+nIsh(i)+nSsh(i)+nDel(i)
  ntri = ntri+nBas(i)*(nBas(i)+1)/2
  nSQ = nSQ+nBas(i)**2
  nVV = nVV+nSsh(i)**2
  nBmx = max(nBmx,nBas(i))
end do
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

NCMO = nSQ
call mma_allocate(CMO,NCMO,2,label='LCMO')
CMO(:,1) = CMOI(1:NCMO)

nOA = 0
do iSym=1,nSym  ! setup info
  lnFro(iSym) = nFro(iSym)
  lnOcc(iSym) = nIsh(iSym)
  nOA = nOA+lnOcc(iSym)
  lnVir(iSym) = nSsh(iSym)
  lnOrb(iSym) = lnOcc(iSym)+lnVir(iSym)
  lnDel(iSym) = nDel(iSym)
end do

call mma_allocate(OrbE,nOrb,4,label='Eorb')
jOff = 1
kOff = 1
lOff = 1
do iSym=1,nSym
  jp = lOff+nFro(iSym)
  jOcc = jOff
  call dcopy_(nIsh(iSym),EOcc(jOcc),1,OrbE(jp,1),1)
  jVir = kOff
  jp = jp+nIsh(iSym)
  call dcopy_(nSsh(iSym),EVir(jVir),1,OrbE(jp,1),1)
  jOff = jOff+nIsh(iSym)
  kOff = kOff+nSsh(iSym)
  lOff = lOff+nBas(iSym)
end do
iOff = 1
jOff = 1
kOff = 1
do iSym=1,nSym
  ifr = iOff+nFro(iSym)
  ito = jOff
  call dcopy_(nIsh(iSym),OrbE(ifr,1),1,OrbE(ito,3),1)
  ifr = iOff+nFro(iSym)+nIsh(iSym)
  ito = kOff
  call dcopy_(nSsh(iSym),OrbE(ifr,1),1,OrbE(ito,4),1)
  iOff = iOff+nBas(iSym)
  jOff = jOff+nIsh(iSym)
  kOff = kOff+nSsh(iSym)
end do
call mma_allocate(X,nVV+nOA,label='Dmat')
X(:) = Zero

call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
CMO(:,2) = Zero
iOff = 1
do iSym=1,nSym
  kfr = iOff+nBas(iSym)*nFro(iSym)
  kto = iOff+nBas(iSym)*lnFro(iSym)
  call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr,1),1,CMO(kto,2),1)
  kfr = iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr,1),1,CMO(kto,2),1)
  iOff = iOff+nBas(iSym)**2
end do
call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,CMO(:,2),OrbE(:,3),OrbE(:,4),X(1:nVV),X(nVV+1:))
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

! Diagonalize the pseudodensity to get natural virtual orbitals
! -------------------------------------------------------------
iOff = 1
jOff = 1
do iSym=1,nSym
  if (nSsh(iSym) > 0) then
    jD = iOff
    ! Eigenvectors will be in increasing order of eigenvalues
    call Eigen_Molcas(nSsh(iSym),X(jD),OrbE(:,2),OrbE(:,1))
    ! Reorder to get relevant eigenpairs first
    do j=1,nSsh(iSym)/2
      do i=1,nSsh(iSym)
        lij = jD-1+nSsh(iSym)*(j-1)+i
        kij = jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
        tmp = X(lij)
        X(lij) = X(kij)
        X(kij) = tmp
      end do
      tmp = OrbE(j,2)
      OrbE(j,2) = OrbE(nSsh(iSym)-j+1,2)
      OrbE(nSsh(iSym)-j+1,2) = tmp
    end do

    ! Compute new MO coeff. : X=C*U
    kfr = jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    kto = jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,CMO(kfr,2),nBas(iSym),X(jD),nSsh(iSym),Zero,CMO(kto,1),nBas(iSym))
    iOff = iOff+nSsh(iSym)**2
    TrDF(iSym) = ddot_(nSsh(iSym),OrbE(:,2),1,[One],0)
    ns_V(iSym) = int(vfrac*nSsh(iSym))
    TrDP(iSym) = ddot_(ns_V(iSym),OrbE(:,2),1,[One],0)
  end if
  jOff = jOff+nBas(iSym)**2
end do
write(u6,*) '------------------------------------------------------'
write(u6,*) '   Symm.     Trace     (Full Dmat)     (Partial Dmat) '
write(u6,*) '------------------------------------------------------'
STrDF = Zero
STrDP = Zero
do iSym=1,nSym
  write(u6,'(4X,I4,14X,G13.6,5X,G13.6)') iSym,TrDF(iSym),TrDP(iSym)
  STrDF = STrDF+TrDF(iSym)
  STrDP = STrDP+TrDP(iSym)
end do
write(u6,*) '------------------------------------------------------'
write(u6,'(A,G13.6,5X,G13.6)') '   Sum :              ',STrDF,STrDP
write(u6,*) '------------------------------------------------------'

! Update the nSsh, nDel for FNO-MP2
nSsh_t = 0
do iSym=1,nSym
  lnOrb(iSym) = lnOrb(iSym)-nSsh(iSym)+ns_V(iSym)
  nDel(iSym) = nDel(iSym)+nSsh(iSym)-ns_V(iSym)
  nSsh(iSym) = ns_V(iSym)
  nSsh_t = nSsh_t+nSsh(iSym)
  nAuxO(iSym) = nBas(iSym)-nDel(iSym)
end do
call Put_iArray('nDelPT',nDel,nSym)
call Put_iArray('nOrb',nAuxO,nSym)

! Reset MP2_small for this new call to ChoMP2_Drv
call Check_Amp2(nSym,lnOcc,nSsh,iSkip)
MP2_small = iSkip > 0
if (MP2_small) then

  call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh)

  call mma_allocate(iD,nOrb,label='iD_orb')
  do k=1,nOrb
    iD(k) = k
  end do
  lOff = 1
  kOff = 1
  jOff = 1
  iOff = 1
  do iSym=1,nSym  ! canonical orb. in the reduced virtual space
    jD = iOff
    call Get_Can_Lorb(OrbE(lOff,4),OrbE(jOff,1),nSsh(iSym),lnVir(iSym),iD,X(jD))

    kfr = kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    kto = kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    nBx = max(1,nBas(iSym))
    nSx = max(1,nSsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,CMO(kfr,1),nBx,X(jD),nSx,Zero,CMOI(kto),nBx)

    lOff = lOff+lnVir(iSym)
    kOff = kOff+nBas(iSym)**2
    jOff = jOff+nSsh(iSym)
    iOff = iOff+lnVir(iSym)**2
  end do
  call mma_deallocate(iD)

  ! Copy the new Evir to output array
  call dcopy_(nSsh_t,OrbE(:,1),1,EVir,1)

  write(u6,*)
  write(u6,'(A,8I4)') ' Secondary orbitals after selection:',(nSsh(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Deleted orbitals after selection:  ',(nDel(i),i=1,nSym)
  write(u6,*)
  write(u6,*) 'Energies of the active virtual orbitals '
  ii = 0
  do iSym=1,nSym
    if (nSsh(iSym) /= 0) then
      write(u6,*)
      write(u6,'(A,I2,(T40,5F14.6))') ' symmetry species',iSym,(EVir(ii+k),k=1,nSsh(iSym))
      ii = ii+nSsh(iSym)
    end if
  end do
  write(u6,*)

  EMP2 = DeMP2
  DeMP2 = Zero
  if (DoMP2) call ChoMP2_Drv(irc,Dummy,CMOI,OrbE(:,3),OrbE(:,1),X(1:nVV),X(nVV+1:))
  if (irc /= 0) then
    write(u6,*) 'MP2 in truncated virtual space failed !'
    call Abend()
  end if
  EMP2 = -One*(EMP2-DeMP2)
  if (DoMP2) then
    DeMP2 = EMP2
    !write(u6,*)
    !write(u6,'(1x,a,f18.10,a)')'FNO correction:       ',EMP2,'   (estimate)   '
    !write(u6,*)
  else
    EMP2 = Zero
  end if
else
  write(u6,*)
  write(u6,*) 'We found ZERO amplitudes T(ai,bj) with the final '
  write(u6,*) 'combinations of inactive and virtual orbitals !! '
  write(u6,*) 'Check your input and rerun the calculation! Bye!!'
  call Abend()
end if

! Update runfile for subsequent calcs (e.g., CHCC)
iOff = 1
jOff = 1
kOff = 1
do iSym=1,nSym
  ifr = iOff
  ito = kOff+nFro(iSym)
  call dcopy_(nIsh(iSym),EOcc(ifr),1,OrbE(ito,1),1)
  ifr = jOff
  ito = ito+nIsh(iSym)
  call dcopy_(nSsh(iSym),EVir(ifr),1,OrbE(ito,1),1)
  iOff = iOff+nIsh(iSym)
  jOff = jOff+nSsh(iSym)
  kOff = kOff+nBas(iSym)
end do
call Put_dArray('OrbE',OrbE(:,1),nOrb)
call Put_dArray('Last orbitals',CMOI,NCMO)

call mma_deallocate(X)
call mma_deallocate(OrbE)

call mma_deallocate(CMO)

return

end subroutine FNO_MP2
