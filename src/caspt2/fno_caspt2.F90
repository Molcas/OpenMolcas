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

subroutine FNO_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,vfrac,IFQCAN,DoMP2,EMP2,CMO,NCMO)
!***********************************************************************
!                                                                      *
!     Purpose:  setup of Frozen Natural Orbitals CASPT2 (FNO-CASPT2).  *
!                                                                      *
!     Author:   F. Aquilante  (Geneva, May  2008)                      *
!                                                                      *
!***********************************************************************

use InputData, only: Input
use ChoMP2, only: DeMP2, MP2_small, shf
use Molcas, only: MxBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), NCMO
integer(kind=iwp), intent(inout) :: nSsh(nSym), nDel(nSym), IFQCAN
real(kind=wp), intent(in) :: vfrac
logical(kind=iwp), intent(in) :: DoMP2
real(kind=wp), intent(inout) :: EMP2, CMO(NCMO)
integer(kind=iwp) :: i, iAoff, ifr, ioff, ip_X, ip_Y, ipOrbE, iSkip, iSym, ito, j, jD, joff, k, kfr, kij, koff, kto, lij, &
                     lnDel(8), lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), lOff, nAct(8), nBasT, nBx, nOA, nOrb, ns_V(8), nSx, nVV
real(kind=wp) :: Delta_TrD, Dummy, STrDF, STrDP, tmp, TrDF(8), TrDP(8)
integer(kind=iwp), allocatable :: ID(:)
real(kind=wp), allocatable :: CMOX(:,:), DMAT(:), OrbE(:,:)

irc = 0
MP2_small = .false.
shf = Input%RegFNO

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nAct(:) = 0
ns_V(:) = 0
nBasT = sum(nBas(:))
nOrb = sum(nFro(:)+nIsh(:)+nAsh(:)+nSsh(:)+nDel(:))
nVV = sum(nSsh(:)**2)
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

!----------------------------------------------------------------------*
!     Read the molecular orbitals from JobIph                          *
!----------------------------------------------------------------------*
if (IFQCAN == 0) write(u6,'(/6X,A)') 'No pseudocanonical RASSCF orbitals found! I will proceed with FDIAG values.'
call mma_allocate(CMOX,NCMO,2,Label='CMOX')
! This is not the best solution, but I wanted to avoid having to rewrite
! the indexing code below just to use the CMO array directly
CMOX(:,1) = CMO(:)
call mma_allocate(OrbE,nOrb,4,Label='OrbE')
call Get_darray('RASSCF OrbE',OrbE(:,1),nOrb)
iAoff = 0
do iSym=1,nSym
  ipOrbE = iAoff+nFro(iSym)+nIsh(iSym)
  nAct(iSym) = count(OrbE(ipOrbE+1:ipOrbE+nAsh(iSym),1) < Zero)
  iAoff = iAoff+nBas(iSym)
end do

! setup info
lnOrb(1:nSym) = nBas(:)
lnFro(1:nSym) = nFro(:)
lnOcc(1:nSym) = nIsh(:)+nAct(1:nSym)
lnVir(1:nSym) = nSsh(:)
lnDel(1:nSym) = nDel(:)
nOA = sum(lnOcc(1:nSym))

ioff = 1
joff = 1
koff = 1
do iSym=1,nSym
  ifr = ioff+nFro(iSym)
  ito = joff
  OrbE(ito:ito+lnOcc(iSym)-1,3) = OrbE(ifr:ifr+lnOcc(iSym)-1,1)
  ifr = ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
  ito = koff
  OrbE(ito:ito+nSsh(iSym)-1,4) = OrbE(ifr:ifr+nSsh(iSym)-1,1)
  ioff = ioff+nBas(iSym)
  joff = joff+lnOcc(iSym)
  koff = koff+nSsh(iSym)
end do
call mma_allocate(Dmat,nVV+nOA,LABEL='DMAT')
DMAT(:) = Zero
ip_X = 1
ip_Y = ip_X+nVV

call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
CMOX(:,2) = Zero
iOff = 0
do iSym=1,nSym
  kfr = 1+iOff+nBas(iSym)*nFro(iSym)
  kto = 1+iOff+nBas(iSym)*lnFro(iSym)
  CMOX(kto:kto+nBas(iSym)*lnOcc(iSym)-1,2) = CMOX(kfr:kfr+nBas(iSym)*lnOcc(iSym)-1,1)
  kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  CMOX(kto:kto+nBas(iSym)*lnVir(iSym)-1,2) = CMOX(kfr:kfr+nBas(iSym)*lnVir(iSym)-1,1)
  iOff = iOff+nBas(iSym)**2
end do
call Check_Amp(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,CMOX(:,2),OrbE(:,3),OrbE(:,4),DMAT(ip_X),DMAT(ip_Y))
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
iOff = 0
jOff = 0
do iSym=1,nSym
  if (nSsh(iSym) > 0) then
    jD = ip_X+iOff
    ! Eigenvectors will be in increasing order of eigenvalues
    call Eigen_Molcas(nSsh(iSym),DMAT(jD),OrbE(:,2),OrbE(:,1))
    ! Reorder to get relevant eigenpairs first
    do j=1,nSsh(iSym)/2
      do i=1,nSsh(iSym)
        lij = jD-1+nSsh(iSym)*(j-1)+i
        kij = jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
        tmp = DMAT(lij)
        DMAT(lij) = DMAT(kij)
        DMAT(kij) = tmp
      end do
      tmp = OrbE(j,2)
      OrbE(j,2) = OrbE(nSsh(iSym)-j+1,2)
      OrbE(nSsh(iSym)-j+1,2) = tmp
    end do

    ! Compute new MO coeff. : X=C*U
    kfr = 1+jOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
    kto = 1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,CMOX(kfr,2),nBas(iSym),DMAT(jD),nSsh(iSym),Zero,CMOX(kto,1),nBas(iSym))
    iOff = iOff+nSsh(iSym)**2
    TrDF(iSym) = sum(OrbE(1:nSsh(iSym),2))
    if (vfrac >= Zero) then
      ns_V(iSym) = int(vfrac*nSsh(iSym))
      TrDP(iSym) = sum(OrbE(1:ns_V(iSym),2))
    else
      ns_V(iSym) = nSsh(iSym)-1
      TrDP(iSym) = sum(OrbE(1:ns_V(iSym),2))
      Delta_TrD = TrDP(iSym)-TrDF(iSym) ! this is negative
      Delta_TrD = Delta_TrD/TrDF(iSym)
      do while (Delta_TrD > vfrac)
        ns_V(iSym) = ns_V(iSym)-1
        TrDP(iSym) = sum(OrbE(1:ns_V(iSym),2))
        Delta_TrD = (TrDP(iSym)-TrDF(iSym))/TrDF(iSym)
      end do
    end if
  end if
  jOff = jOff+nBas(iSym)**2
end do
write(u6,*) '------------------------------------------------------'
write(u6,*) '   Symm.     Trace     (Full Dmat)     (Partial Dmat) '
write(u6,*) '------------------------------------------------------'
STrDF = Zero
STrDP = Zero
do iSym=1,nSym
  write(u6,'(4X,I4,15X,G13.6,4X,G13.6)') iSym,TrDF(iSym),TrDP(iSym)
  STrDF = STrDF+TrDF(iSym)
  STrDP = STrDP+TrDP(iSym)
end do
write(u6,*) '------------------------------------------------------'
write(u6,'(A,G13.6,4X,G13.6)') '   Sum :               ',STrDF,STrDP
write(u6,*) '------------------------------------------------------'

! Update the nSsh, nDel for FNO-CASPT2
nDel(:) = nDel(:)+nSsh(:)-ns_V(1:nSym)
nSsh(:) = ns_V(1:nSym)

! Write the resorted MOs back to JobIph

IFQCAN = 0 ! MOs need to be recanonicalized on exit
CMO(:) = CMOX(:,1)

! Reset MP2_small for this new call to ChoMP2_Drv
call Check_Amp(nSym,lnOcc,nSsh,iSkip)
MP2_small = DoMP2 .and. (iSkip > 0)
if (MP2_small) then

  call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh)

  call mma_allocate(iD,nOrb,Label='iD')
  iD(:) = [(k,k=1,nOrb)]
  lOff = 1
  kOff = 1
  jOff = 1
  iOff = 0
  do iSym=1,nSym  ! canonical orb. in the reduced virtual space
    jD = ip_X+iOff
    call Get_Can_Lorb(OrbE(lOff,4),OrbE(jOff,1),nSsh(iSym),lnVir(iSym),iD,DMAT(jD))

    kfr = kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    kto = kOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
    nBx = max(1,nBas(iSym))
    nSx = max(1,nSsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,CMOX(kfr,1),nBx,DMAT(jD),nSx,Zero,CMOX(kto,2),nBx)

    lOff = lOff+lnVir(iSym)
    kOff = kOff+nBas(iSym)**2
    jOff = jOff+nSsh(iSym)
    iOff = iOff+lnVir(iSym)**2
  end do
  call mma_deallocate(iD)

  EMP2 = DeMP2
  DeMP2 = Zero
  call ChoMP2_Drv(irc,Dummy,CMOX(:,2),OrbE(:,3),OrbE(:,1),DMAT(ip_X),DMAT(ip_Y))
  if (irc /= 0) then
    write(u6,*) 'MP2 in truncated virtual space failed !'
    call Abend()
  end if
  EMP2 = -One*(EMP2-DeMP2)
end if
call mma_deallocate(Dmat)
call mma_deallocate(OrbE)

call mma_deallocate(CMOX)

end subroutine FNO_CASPT2
