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
integer(kind=iwp) :: i, iCMO, ifr, ii, ioff, ip_iD, ip_X, ip_Y, ip_Z, ip_ZZ, ipEorb, ipOrbE, iSkip, iSym, ito, j, jD, jOcc, jOff, &
                     jp, jVir, k, kEOcc, kEVir, kfr, kij, kOff, kto, LCMO, lij, lnDel(8), lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), &
                     lOff, nAuxO(8), nBasT, nBmx, nBx, NCMO, nOA, nOrb, ns_V(8), nSQ, nSsh_t, nSx, ntri, nVV
real(kind=wp) :: Dummy, STrDF, STrDP, tmp, TrDF(8), TrDP(8)
real(kind=wp), external :: ddot_
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "chfnopt.fh"

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
call GETMEM('LCMO','ALLO','REAL',LCMO,2*NCMO)
iCMO = LCMO+NCMO
call DCOPY_(NCMO,CMOI,1,WORK(LCMO),1)

nOA = 0
do iSym=1,nSym  ! setup info
  lnFro(iSym) = nFro(iSym)
  lnOcc(iSym) = nIsh(iSym)
  nOA = nOA+lnOcc(iSym)
  lnVir(iSym) = nSsh(iSym)
  lnOrb(iSym) = lnOcc(iSym)+lnVir(iSym)
  lnDel(iSym) = nDel(iSym)
end do

call GetMem('Eorb','Allo','Real',ipOrbE,4*nOrb)
jOff = 0
kOff = 0
lOff = 0
do iSym=1,nSym
  jp = ipOrbE+lOff+nFro(iSym)
  jOcc = jOff+1
  call dcopy_(nIsh(iSym),EOcc(jOcc),1,Work(jp),1)
  jVir = kOff+1
  jp = jp+nIsh(iSym)
  call dcopy_(nSsh(iSym),EVir(jVir),1,Work(jp),1)
  jOff = jOff+nIsh(iSym)
  kOff = kOff+nSsh(iSym)
  lOff = lOff+nBas(iSym)
end do
ip_ZZ = ipOrbE
ipEorb = ipOrbE+nOrb
ip_Z = ipEorb
kEOcc = ipEorb+nOrb
kEVir = kEOcc+nOrb
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = ipOrbE+ioff+nFro(iSym)
  ito = kEOcc+joff
  call dcopy_(nIsh(iSym),Work(ifr),1,Work(ito),1)
  ifr = ipOrbE+ioff+nFro(iSym)+nIsh(iSym)
  ito = kEVir+koff
  call dcopy_(nSsh(iSym),Work(ifr),1,Work(ito),1)
  ioff = ioff+nBas(iSym)
  joff = joff+nIsh(iSym)
  koff = koff+nSsh(iSym)
end do
call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
ip_Y = ip_X+nVV
call FZero(Work(ip_X),nVV+nOA)

call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y)
call FZero(Work(iCMO),NCMO)
iOff = 0
do iSym=1,nSym
  kfr = LCMO+iOff+nBas(iSym)*nFro(iSym)
  kto = iCMO+iOff+nBas(iSym)*lnFro(iSym)
  call dcopy_(nBas(iSym)*lnOcc(iSym),Work(kfr),1,Work(kto),1)
  kfr = LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  call dcopy_(nBas(iSym)*lnVir(iSym),Work(kfr),1,Work(kto),1)
  iOff = iOff+nBas(iSym)**2
end do
call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
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
    call Eigen_Molcas(nSsh(iSym),Work(jD),Work(ip_Z),Work(ip_ZZ))
    ! Reorder to get relevant eigenpairs first
    do j=1,nSsh(iSym)/2
      do i=1,nSsh(iSym)
        lij = jD-1+nSsh(iSym)*(j-1)+i
        kij = jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
        tmp = Work(lij)
        Work(lij) = Work(kij)
        Work(kij) = tmp
      end do
      tmp = Work(ip_Z-1+j)
      Work(ip_Z-1+j) = Work(ip_Z+nSsh(iSym)-j)
      Work(ip_Z+nSsh(iSym)-j) = tmp
    end do

    ! Compute new MO coeff. : X=C*U
    kfr = iCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    kto = LCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,Work(kfr),nBas(iSym),Work(jD),nSsh(iSym),Zero,Work(kto),nBas(iSym))
    iOff = iOff+nSsh(iSym)**2
    TrDF(iSym) = ddot_(nSsh(iSym),Work(ip_Z),1,[One],0)
    ns_V(iSym) = int(vfrac*nSsh(iSym))
    TrDP(iSym) = ddot_(ns_V(iSym),Work(ip_Z),1,[One],0)
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

  call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh,ip_X,ip_Y)

  call GetMem('iD_orb','Allo','Inte',ip_iD,nOrb)
  do k=1,nOrb
    iWork(ip_iD-1+k) = k
  end do
  lOff = 0
  kOff = 0
  jOff = 0
  iOff = 0
  do iSym=1,nSym  ! canonical orb. in the reduced virtual space
    jD = ip_X+iOff
    call Get_Can_Lorb(Work(kEVir+lOff),Work(ipOrbE+jOff),nSsh(iSym),lnVir(iSym),iWork(ip_iD),Work(jD),iSym)

    kfr = LCMO+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    kto = 1+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
    nBx = max(1,nBas(iSym))
    nSx = max(1,nSsh(iSym))
    call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),One,Work(kfr),nBx,Work(jD),nSx,Zero,CMOI(kto),nBx)

    lOff = lOff+lnVir(iSym)
    kOff = kOff+nBas(iSym)**2
    jOff = jOff+nSsh(iSym)
    iOff = iOff+lnVir(iSym)**2
  end do
  call GetMem('iD_orb','Free','Inte',ip_iD,nOrb)
  kEVir = ipOrbE

  ! Copy the new Evir to output array
  call dcopy_(nSsh_t,Work(kEVir),1,EVir,1)

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
  if (DoMP2) call ChoMP2_Drv(irc,Dummy,CMOI,Work(kEOcc),Work(kEVir))
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
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = 1+ioff
  ito = ipOrbE+koff+nFro(iSym)
  call dcopy_(nIsh(iSym),EOcc(ifr),1,Work(ito),1)
  ifr = 1+joff
  ito = ito+nIsh(iSym)
  call dcopy_(nSsh(iSym),EVir(ifr),1,Work(ito),1)
  ioff = ioff+nIsh(iSym)
  joff = joff+nSsh(iSym)
  koff = koff+nBas(iSym)
end do
call Put_dArray('OrbE',Work(ipOrbE),nOrb)
call Put_dArray('Last orbitals',CMOI,NCMO)

call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
call GetMem('Eorb','Free','Real',ipOrbE,4*nOrb)

call GETMEM('LCMO','FREE','REAL',LCMO,2*NCMO)

return

end subroutine FNO_MP2
