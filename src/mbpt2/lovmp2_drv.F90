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

subroutine LovMP2_Drv(irc,EMP2,CMO,EOcc,EVir,NamAct,n_Acta,Thrs,Do_MP2,allVir)
! Purpose: setup of Localized occupied-virtual MP2  (LovMP2)
!          The MP2 correction to the energy will later be computed
!          only for the "active region" of the molecule.
!          The MP2 correction due to the remaining frozen region
!          is computed here if Do_MP2=.true.
!
! Author:  F. Aquilante  (Geneva, Jun. 2008)

use MBPT2_Global, only: nBas
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(inout) :: CMO(*), EOcc(*), EVir(*)
integer(kind=iwp), intent(in) :: n_Acta
character(len=LenIn), intent(in) :: NamAct(n_Acta)
real(kind=wp), intent(in) :: Thrs
logical(kind=iwp), intent(in) :: Do_MP2, allVir
integer(kind=iwp) :: i, ia, iComp, iDo, ie, ifr, ii, ik, iloc, iOff, iOpt, iSkip, iSym, isymlbl, ito, iV, ja, jDo, jk, jloc, jOff, &
                     k, ka, kfr, kk, kOff, kto, lnDel(8), lnDel2(8), lnFro(8), lnFro2(8), lnOrb(8), lnOcc(8), lnOcc2(8), lnVir(8), &
                     lnVir2(8), lOff, lsq, ltri, nAuxO(8), nBmx, nOA, ns_O(8), ns_V(8), nSQ, ntri, nVV, nxBasT, nxOrb, nZero(8)
real(kind=wp) :: Dummy, EFRO, EOSF, StrA, STrF, STrX, Thrd, TrA(8), TrF(8), TrX(8)
logical(kind=iwp) :: ortho
character(len=8) :: Label
integer(kind=iwp), allocatable :: iD_vir(:)
real(kind=wp), allocatable :: EOrb(:,:), LCMO(:,:), S(:), Saa(:), SQ(:), X(:)
character(len=LenIn8), allocatable :: UBName(:)
real(kind=wp), external :: ddot_
#include "corbinf.fh"
#include "chomp2_cfg.fh"

irc = 0
EMP2 = Zero
EFRO = Zero
EOSF = Zero
iDo = 0
jDo = 0

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nxBasT = 0
ntri = 0
nSQ = 0
nBmx = 0
nxOrb = 0
do i=1,nSym
  TrA(i) = Zero
  TrF(i) = Zero
  TrX(i) = Zero
  nZero(i) = 0
  nxBasT = nxBasT+nBas(i)
  nxOrb = nxOrb+nFro(i)+nOcc(i)+nExt(i)+nDel(i)
  ntri = ntri+nBas(i)*(nBas(i)+1)/2
  nSQ = nSQ+nBas(i)**2
  nBmx = max(nBmx,nBas(i))
end do
if (nxBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

call mma_allocate(UBName,nxBasT,label='UBName')
call Get_cArray('Unique Basis Names',UBName,(LenIn8)*nxBasT)

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
call mma_allocate(SQ,nSQ,label='SMAT')
call mma_allocate(S,nTri,label='SLT')
isymlbl = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,S,isymlbl)
if (irc /= 0) then
  return
end if
ltri = 1
lsq = 1
do iSym=1,nSym
  call Square(S(ltri),SQ(lsq),1,nBas(iSym),nBas(iSym))
  ltri = ltri+nBas(iSym)*(nBas(iSym)+1)/2
  lsq = lsq+nBas(iSym)**2
end do
call mma_deallocate(S)

write(u6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
write(u6,*)
if (n_Acta /= 0) then
  write(u6,'(A,I3,A)') ' Selected ',n_Acta,' atoms: '
  write(u6,*)
  write(u6,*) (NamAct(i),i=1,n_Acta)
  write(u6,*)
else if (.not. Do_MP2) then
  write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
  call finalize()
  return
else
  write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
end if

!----------------------------------------------------------------------*
call Get_Tr_Dab(nSym,nBas,nFro,nOcc,nExt,nDel,CMO,EOcc,EVir,TrX)
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!     Localize the inactive and virtual orbitals                       *
!                                                                      *
!        1) inactive orbitals ---> cholesky orbitals (orthonormal)     *
!        2) virtual orbitals ---> lin. indep. PAOs (non-orthonormal)   *
!----------------------------------------------------------------------*
call mma_allocate(LCMO,nSQ,3,label='LCMO')
LCMO(:,2) = CMO(1:nSQ)
Thrd = 1.0e-6_wp
call mma_allocate(iD_vir,nxBasT,label='ID_vir')
call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nOcc,nZero,nExt,LCMO(:,2),SQ,iD_vir)

if (irc /= 0) then
  write(u6,*) 'Localization failed in LovMP2'
  call Abend()
end if

call mma_allocate(EOrb,nxOrb,2,label='Eorb')
EOrb(:,1) = EOcc(1:nxOrb)
EOrb(:,2) = EVir(1:nxOrb)

call mma_allocate(Saa,nxOrb,label='Saa')
Saa(:) = One

!----------------------------------------------------------------------*
!     Occupied orbital selection                                       *
!----------------------------------------------------------------------*
iOff = 1
kOff = 1
do iSym=1,nSym
  jOff = iOff+nBas(iSym)*nFro(iSym)
  call dcopy_(nBas(iSym)*nOcc(iSym),LCMO(jOff,2),1,LCMO(kOff,3),1)
  call dcopy_(nBas(iSym)*nOcc(iSym),CMO(jOff),1,LCMO(kOff,1),1)
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nOcc(iSym)
end do
ortho = .true.

call get_Orb_select(irc,LCMO(:,1),LCMO(:,3),EOrb(:,1),SQ,Saa,UBName,NamAct,nSym,n_Acta,nOcc,nBas,ortho,Thrs,ns_O)
if (irc /= 0) then
  return
end if
iOff = 1
kOff = 1
do iSym=1,nSym
  lOff = iOff+nBas(iSym)*nFro(iSym)
  do ik=nOcc(iSym),1,-1
    jOff = kOff+nBas(iSym)*(ik-1)
    call dcopy_(nBas(iSym),LCMO(jOff,1),1,CMO(lOff),1)
    lOff = lOff+nBas(iSym)
  end do
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nOcc(iSym)
end do
iloc = 0
loff = 0
do iSym=1,nSym
  do ik=nOcc(iSym),1,-1
    ie = loff+ik
    iloc = iloc+1
    EOcc(iloc) = EOrb(ie,1)
  end do
  loff = loff+nOcc(iSym)
end do

if (allVir) then
  do iSym=1,nSym
    ns_V(iSym) = nExt(iSym)
  end do
else
  !--------------------------------------------------------------------*
  !     Virtual orbital selection                                      *
  !--------------------------------------------------------------------*
  iOff = 1
  kOff = 1
  do iSym=1,nSym
    jOff = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
    call dcopy_(nBas(iSym)*nExt(iSym),LCMO(jOff,2),1,LCMO(kOff,3),1)
    call dcopy_(nBas(iSym)*nExt(iSym),CMO(jOff),1,LCMO(kOff,1),1)
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*nExt(iSym)
  end do
  ortho = .false.

  call get_Vir_select(irc,LCMO(:,1),LCMO(:,3),EOrb(:,2),SQ,UBName,NamAct,iD_vir,nSym,n_Acta,nExt,nBas,ortho,ns_V)
  if (irc /= 0) then
    return
  end if
  call mma_deallocate(iD_vir)
  iOff = 1
  kOff = 1
  do iSym=1,nSym
    jOff = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
    call dcopy_(nBas(iSym)*nExt(iSym),LCMO(kOff,1),1,CMO(jOff),1)
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*nExt(iSym)
  end do
  iloc = 0
  loff = 0
  do iSym=1,nSym
    do ik=1,nExt(iSym)
      ie = loff+ik
      iloc = iloc+1
      EVir(iloc) = EOrb(ie,2)
    end do
    loff = loff+nExt(iSym)
  end do
end if

call mma_deallocate(UBName)

iDo = 0
jDo = 0
do iSym=1,nSym  ! setup info
  lnOcc2(iSym) = nOcc(iSym)
  lnFro2(iSym) = nFro(iSym)
  lnDel2(iSym) = nDel(iSym)
  lnVir2(iSym) = nExt(iSym)
  lnOrb(iSym) = nBas(iSym)
  lnOcc(iSym) = nOcc(iSym)-ns_O(iSym)
  lnFro(iSym) = nFro(iSym)+ns_O(iSym)
  lnDel(iSym) = nDel(iSym)+ns_V(iSym)
  lnVir(iSym) = nExt(iSym)-ns_V(iSym)
  iDo = max(iDo,lnOcc(iSym))
  jDo = max(jDo,lnVir(iSym))
end do
if (min(iDo,jDo) /= 0) then
  !--------------------------------------------------------------------*
  !     MP2 calculation on the Frozen region                           *
  !--------------------------------------------------------------------*
  if (Do_MP2) then

    iloc = 1
    jloc = 1
    loff = 0
    joff = 0
    do iSym=1,nSym
      do ik=ns_V(iSym)+1,nExt(iSym)
        ie = loff+ik
        EOrb(iloc,2) = EVir(ie)
        iloc = iloc+1
      end do
      loff = loff+nExt(iSym)
      do ik=1,lnOcc(iSym)
        ie = joff+ik
        EOrb(jloc,1) = EOcc(ie)
        jloc = jloc+1
      end do
      joff = joff+nOcc(iSym)
    end do
    iOff = 1
    nVV = 0
    nOA = 0
    do iSym=1,nSym
      kfr = iOff+nBas(iSym)*nFro(iSym)
      kto = iOff+nBas(iSym)*lnFro(iSym)
      call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,LCMO(kto,1),1)
      kfr = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym)+ns_V(iSym))
      kto = kto+nBas(iSym)*lnOcc(iSym)
      call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,LCMO(kto,1),1)
      iOff = iOff+nBas(iSym)**2
      nVV = nVV+lnVir(iSym)**2
      nOA = nOA+lnOcc(iSym)
    end do
    call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
    if (iSkip > 0) then
      call mma_allocate(X,nVV+nOA,label='Dmat')
      X(:) = Zero
      call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,X(1:nVV),X(nVV+1:),.true.)
      call ChoMP2_Drv(irc,Dummy,LCMO(:,1),EOrb(:,1),EOrb(:,2))
      call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,X(1:nVV),X(nVV+1:),.false.) ! compute energy and not Dab
      call ChoMP2_Drv(irc,EFRO,LCMO(:,1),EOrb(:,1),EOrb(:,2))
      if (irc /= 0) then
        write(u6,*) 'Frozen region MP2 failed'
        call Abend()
      else
        write(u6,'(A,F20.10,A)') ' Frozen region E2 contrib. = ',EFRO,' a.u.'
        EOSF = EOSMP2
        write(u6,'(A,F20.10,A)') ' (Opposite-Spin contrib.   = ',-EOSF,' )'
        write(u6,*)
      end if
      iV = 1
      do iSym=1,nSym
        TrF(iSym) = ddot_(lnVir(iSym),X(iV),1+lnVir(iSym),[One],0)
        iV = iV+lnVir(iSym)**2
      end do
      call mma_deallocate(X)
      do iSym=1,nSym
        nOcc(iSym) = lnOcc2(iSym)
        nFro(iSym) = lnFro2(iSym)
        nDel(iSym) = lnDel2(iSym)
        nExt(iSym) = lnVir2(iSym)
      end do
    end if

  end if
end if

! Update the nFro, nOcc, nExt, nDel for the Active site MP2
do iSym=1,nSym
  nFro(iSym) = nFro(iSym)+nOcc(iSym)-ns_O(iSym)
  nOcc(iSym) = ns_O(iSym)
  nDel(iSym) = nDel(iSym)+nExt(iSym)-ns_V(iSym)
  nExt(iSym) = ns_V(iSym)
  iDo = max(iDo,nOcc(iSym))
  jDo = max(jDo,nExt(iSym))
  nAuxO(iSym) = nBas(iSym)-nDel(iSym)
end do
call Put_iArray('nFroPT',nFro,nSym)
call Put_iArray('nDelPT',nDel,nSym)
call Put_iArray('nOrb',nAuxO,nSym)

call Check_Amp2(nSym,nOcc,nExt,iSkip)
if (iSkip > 0) then
  iOff = 0
  lOff = 0
  jk = 1
  ja = 1
  nVV = 0
  nOA = 0
  do iSym=1,nSym
    do ik=1,nOcc(iSym)
      kk = ik+iOff+lnOcc(iSym)
      EOcc(jk) = EOcc(kk)
      jk = jk+1
    end do
    iOff = iOff+lnOcc(iSym)+nOcc(iSym)
    do ia=1,nExt(iSym)
      ka = ia+lOff
      EVir(ja) = EVir(ka)
      ja = ja+1
    end do
    lOff = lOff+lnVir(iSym)+nExt(iSym)
    nVV = nVV+nExt(iSym)**2
    nOA = nOA+nOcc(iSym)
  end do

  write(u6,*)
  write(u6,'(A,8I4)') ' Frozen orbitals after selection :  ',(nFro(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Occupied orbitals after selection: ',(nOcc(i),i=1,nSym)
  write(u6,*)
  write(u6,*) 'Energies of the active occupied orbitals '
  ii = 0
  do iSym=1,nSym
    if (nOcc(iSym) /= 0) then
      write(u6,*)
      write(u6,'(A,I2,(T40,5F14.6))') ' symmetry species',iSym,(EOcc(ii+k),k=1,nOcc(iSym))
      ii = ii+nOcc(iSym)
    end if
  end do
  write(u6,*)

  write(u6,*)
  write(u6,'(A,8I4)') ' Secondary orbitals after selection:',(nExt(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Deleted orbitals after selection:  ',(nDel(i),i=1,nSym)
  write(u6,*)
  write(u6,*) 'Energies of the active virtual orbitals '
  ii = 0
  do iSym=1,nSym
    if (nExt(iSym) /= 0) then
      write(u6,*)
      write(u6,'(A,I2,(T40,5F14.6))') ' symmetry species',iSym,(EVir(ii+k),k=1,nExt(iSym))
      ii = ii+nExt(iSym)
    end if
  end do
  write(u6,*)

  call mma_allocate(X,nVV+nOA,label='Dmat')
  X(:) = Zero
  call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,X(1:nVV),X(nVV+1:),.true.)
  call ChoMP2_Drv(irc,Dummy,CMO,EOcc,EVir)
  call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,X(1:nVV),X(nVV+1:),.false.)
  Wref = Zero
  call ChoMP2_Drv(irc,EMP2,CMO,EOcc,EVir)
  if (irc /= 0) then
    write(u6,*) 'LovMP2 failed'
    call Abend()
  end if
  iV = 1
  do iSym=1,nSym
    TrA(iSym) = ddot_(nExt(iSym),X(iV),1+nExt(iSym),[One],0)
    iV = iV+nExt(iSym)**2
  end do
  call mma_deallocate(X)
end if
write(u6,*) '------------------------------------------------------'
write(u6,*) ' Symm.   Tr(D):  Active        Frozen         Full    '
write(u6,*) '------------------------------------------------------'
STrA = Zero
STrF = Zero
STrX = Zero
do iSym=1,nSym
  write(u6,'(2X,I4,9X,G11.4,3X,G11.4,3X,G11.4)') iSym,TrA(iSym),TrF(iSym),TrX(iSym)
  STrA = STrA+TrA(iSym)
  STrF = STrF+TrF(iSym)
  STrX = STrX+TrX(iSym)
end do
write(u6,*) '------------------------------------------------------'
write(u6,'(A,G11.4,3X,G11.4,3X,G11.4)') '          Sum: ',STrA,STrF,STrX
write(u6,*) '------------------------------------------------------'
write(u6,*)

write(u6,'(A,F20.10,A)') ' Active region E2 contrib. = ',EMP2,' a.u.'
write(u6,'(A,F20.10,A)') ' (Opposite-Spin contrib.   = ',-EOSMP2,' )'
write(u6,*)

XEMP2 = EFRO
EMP2 = EMP2+EFRO
EOSMP2 = EOSMP2+EOSF

! Update runfile for subsequent calcs (e.g., CHCC)
iOff = 1
jOff = 1
kOff = 1
do iSym=1,nSym
  ifr = iOff
  ito = kOff+nFro(iSym)
  call dcopy_(nOcc(iSym),EOcc(ifr),1,EOrb(ito,1),1)
  ifr = jOff
  ito = ito+nOcc(iSym)
  call dcopy_(nExt(iSym),EVir(ifr),1,EOrb(ito,1),1)
  iOff = iOff+nOcc(iSym)
  jOff = jOff+nExt(iSym)
  kOff = kOff+nBas(iSym)
end do
call Put_dArray('OrbE',EOrb(:,1),nxOrb)
call Put_dArray('Last orbitals',CMO,nSQ)

call mma_deallocate(LCMO)
call mma_deallocate(Saa)
call mma_deallocate(EOrb)

call finalize()

return

contains

subroutine finalize()
  if (min(iDo,jDo) == 0) then
    write(u6,*)
    write(u6,*) ' None of the occupied or virtual orbitals has been '
    write(u6,*) ' assigned to the Active region of the molecule.    '
    write(u6,*) ' This is presumably NOT what you want !!!          '
    write(u6,*) ' MP2 will Stop here. Bye Bye !! '
    write(u6,*)
    call Abend()
  end if
  call mma_deallocate(SQ)
end subroutine finalize

end subroutine LovMP2_Drv
