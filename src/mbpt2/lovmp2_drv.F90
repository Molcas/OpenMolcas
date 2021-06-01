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

#include "implicit.fh"
dimension CMO(*), EOcc(*), EVir(*)
real*8 Thrs
logical Do_MP2, allVir
#include "Molcas.fh"
character*(LENIN8) Name(mxBas)
character*(LENIN) NamAct(*)
logical ortho
real*8 TrA(8), TrF(8), TrX(8)
integer ns_O(8), ns_V(8), nZero(8)
integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
integer lnOcc2(8), lnFro2(8), lnDel2(8), lnVir2(8), nAuxO(8)
character*3 ThisNm
character*10 SecNam
parameter(SecNam='LovMP2_Drv',ThisNm='Dry')
#include "itmax.fh"
#include "real.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "WrkSpc.fh"
#include "chomp2_cfg.fh"

irc = 0
EMP2 = 0.0d0
EFRO = 0.0d0
EOSF = 0.0d0
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
  TrA(i) = 0.0d0
  TrF(i) = 0.0d0
  TrX(i) = 0.0d0
  nZero(i) = 0
  nxBasT = nxBasT+nBas(i)
  nxOrb = nxOrb+nFro(i)+nOcc(i)+nExt(i)+nDel(i)
  ntri = ntri+nBas(i)*(nBas(i)+1)/2
  nSQ = nSQ+nBas(i)**2
  nBmx = max(nBmx,nBas(i))
end do
if (nxBasT > mxBas) then
  write(6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

call Get_cArray('Unique Basis Names',Name,(LENIN8)*nxBasT)

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
call GetMem('SMAT','ALLO','REAL',ipSQ,nSQ)
call GetMem('SLT','ALLO','REAL',ipS,nTri)
isymlbl = 1
call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
if (irc /= 0) then
  return
end if
ltri = 0
lsq = 0
do iSym=1,nSym
  call Square(Work(ipS+ltri),Work(ipSQ+lsq),1,nBas(iSym),nBas(iSym))
  ltri = ltri+nBas(iSym)*(nBas(iSym)+1)/2
  lsq = lsq+nBas(iSym)**2
end do
call GetMem('SLT','FREE','REAL',ipS,nTri)

write(6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
write(6,*)
if (n_Acta /= 0) then
  write(6,'(A,I3,A)') ' Selected ',n_Acta,' atoms: '
  write(6,*)
  write(6,*) (NamAct(i),i=1,n_Acta)
  write(6,*)
else if (.not. Do_MP2) then
  write(6,'(A,18A4)') ' Selected atoms: *** None *** '
  Go To 2000
else
  write(6,'(A,18A4)') ' Selected atoms: *** None *** '
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
call GetMem('LCMO','ALLO','REAL',iCMO,3*nSQ)
ipXMO = iCMO+nSQ
iXMO = ipXMO+nSQ
call dcopy_(nSQ,CMO,1,Work(ipXMO),1)
Thrd = 1.0d-06
call GetMem('ID_vir','Allo','Inte',iD_vir,nxBasT)
call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nOcc,nZero,nExt,Work(ipXMO),Work(ipSQ),iWork(iD_vir))

if (irc /= 0) then
  write(6,*) 'Localization failed in LovMP2'
  call Abend()
end if

call GetMem('Eorb','Allo','Real',kEOcc,2*nxOrb)
kEVir = kEOcc+nxOrb
call dcopy_(nxOrb,EOcc,1,Work(kEOcc),1)
call dcopy_(nxOrb,EVir,1,Work(kEVir),1)

call GetMem('Saa','Allo','Real',ipSaa,nxOrb)
call dcopy_(nxOrb,[1.0d0],0,Work(ipSaa),1)

!----------------------------------------------------------------------*
!     Occupied orbital selection                                       *
!----------------------------------------------------------------------*
iOff = 0
kOff = 0
do iSym=1,nSym
  jOff = iOff+nBas(iSym)*nFro(iSym)
  call dcopy_(nBas(iSym)*nOcc(iSym),Work(ipXMO+jOff),1,Work(iXMO+kOff),1)
  call dcopy_(nBas(iSym)*nOcc(iSym),CMO(1+jOff),1,Work(iCMO+kOff),1)
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nOcc(iSym)
end do
ortho = .true.

call get_Orb_select(irc,Work(iCMO),Work(iXMO),Work(kEOcc),Work(ipSQ),Work(ipSaa),Name,NamAct,nSym,n_Acta,nOcc,nBas,ortho,Thrs,ns_O)
if (irc /= 0) then
  return
end if
iOff = 0
kOff = 0
do iSym=1,nSym
  lOff = iOff+nBas(iSym)*nFro(iSym)
  do ik=nOcc(iSym),1,-1
    jOff = kOff+nBas(iSym)*(ik-1)
    call dcopy_(nBas(iSym),Work(iCMO+jOff),1,CMO(1+lOff),1)
    lOff = lOff+nBas(iSym)
  end do
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nOcc(iSym)
end do
iloc = 0
loff = 0
do iSym=1,nSym
  do ik=nOcc(iSym),1,-1
    ie = kEOcc+loff+ik-1
    iloc = iloc+1
    EOcc(iloc) = Work(ie)
  end do
  loff = loff+nOcc(iSym)
end do
!
if (allVir) then
  do iSym=1,nSym
    ns_V(iSym) = nExt(iSym)
  end do
  goto 999
end if
!----------------------------------------------------------------------*
!     Virtual orbital selection                                        *
!----------------------------------------------------------------------*
iOff = 0
kOff = 0
do iSym=1,nSym
  jOff = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
  call dcopy_(nBas(iSym)*nExt(iSym),Work(ipXMO+jOff),1,Work(iXMO+kOff),1)
  call dcopy_(nBas(iSym)*nExt(iSym),CMO(1+jOff),1,Work(iCMO+kOff),1)
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nExt(iSym)
end do
ortho = .false.
!
call get_Vir_select(irc,Work(iCMO),Work(iXMO),Work(kEVir),Work(ipSQ),Name,NamAct,iWork(iD_vir),nSym,n_Acta,nExt,nBas,ortho,ns_V)
if (irc /= 0) then
  return
end if
call GetMem('ID_vir','Free','Inte',iD_vir,nxBasT)
iOff = 0
kOff = 0
do iSym=1,nSym
  jOff = iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
  call dcopy_(nBas(iSym)*nExt(iSym),Work(iCMO+kOff),1,CMO(1+jOff),1)
  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nBas(iSym)*nExt(iSym)
end do
iloc = 0
loff = 0
do iSym=1,nSym
  do ik=1,nExt(iSym)
    ie = kEVir+loff+ik-1
    iloc = iloc+1
    EVir(iloc) = Work(ie)
  end do
  loff = loff+nExt(iSym)
end do

999 continue
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
if (min(iDo,jDo) == 0) goto 1000

!----------------------------------------------------------------------*
!     MP2 calculation on the Frozen region                             *
!----------------------------------------------------------------------*
if (Do_MP2) then

  iloc = 0
  jloc = 0
  loff = 0
  joff = 0
  do iSym=1,nSym
    do ik=ns_V(iSym)+1,nExt(iSym)
      ie = loff+ik
      Work(kEVir+iloc) = EVir(ie)
      iloc = iloc+1
    end do
    loff = loff+nExt(iSym)
    do ik=1,lnOcc(iSym)
      ie = joff+ik
      Work(kEOcc+jloc) = EOcc(ie)
      jloc = jloc+1
    end do
    joff = joff+nOcc(iSym)
  end do
  iOff = 0
  nVV = 0
  nOA = 0
  do iSym=1,nSym
    kfr = 1+iOff+nBas(iSym)*nFro(iSym)
    kto = iCMO+iOff+nBas(iSym)*lnFro(iSym)
    call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
    kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym)+ns_V(iSym))
    kto = kto+nBas(iSym)*lnOcc(iSym)
    call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
    iOff = iOff+nBas(iSym)**2
    nVV = nVV+lnVir(iSym)**2
    nOA = nOA+lnOcc(iSym)
  end do
  call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
  if (iSkip > 0) then
    call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
    ip_Y = ip_X+nVV
    call FZero(Work(ip_X),nVV+nOA)
    call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y,.true.)
    call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
    call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y,.false.) ! compute energy and not Dab
    call ChoMP2_Drv(irc,EFRO,Work(iCMO),Work(kEOcc),Work(kEVir))
    if (irc /= 0) then
      write(6,*) 'Frozen region MP2 failed'
      call Abend()
    else
      write(6,'(A,F20.10,A)') ' Frozen region E2 contrib. = ',EFRO,' a.u.'
      EOSF = EOSMP2
      write(6,'(A,F20.10,A)') ' (Opposite-Spin contrib.   = ',-EOSF,' )'
      write(6,*)
    end if
    iV = ip_X
    do iSym=1,nSym
      TrF(iSym) = ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),[1.0d0],0)
      iV = iV+lnVir(iSym)**2
    end do
    call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
    do iSym=1,nSym
      nOcc(iSym) = lnOcc2(iSym)
      nFro(iSym) = lnFro2(iSym)
      nDel(iSym) = lnDel2(iSym)
      nExt(iSym) = lnVir2(iSym)
    end do
  end if

end if
1000 continue

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

  write(6,*)
  write(6,'(A,8I4)') ' Frozen orbitals after selection :  ',(nFro(i),i=1,nSym)
  write(6,'(A,8I4)') ' Occupied orbitals after selection: ',(nOcc(i),i=1,nSym)
  write(6,*)
  write(6,*) 'Energies of the active occupied orbitals '
  ii = 0
  do iSym=1,nSym
    if (nOcc(iSym) /= 0) then
      write(6,*)
      write(6,'(A,I2,(T40,5F14.6))') ' symmetry species',iSym,(EOcc(ii+k),k=1,nOcc(iSym))
      ii = ii+nOcc(iSym)
    end if
  end do
  write(6,*)

  write(6,*)
  write(6,'(A,8I4)') ' Secondary orbitals after selection:',(nExt(i),i=1,nSym)
  write(6,'(A,8I4)') ' Deleted orbitals after selection:  ',(nDel(i),i=1,nSym)
  write(6,*)
  write(6,*) 'Energies of the active virtual orbitals '
  ii = 0
  do iSym=1,nSym
    if (nExt(iSym) /= 0) then
      write(6,*)
      write(6,'(A,I2,(T40,5F14.6))') ' symmetry species',iSym,(EVir(ii+k),k=1,nExt(iSym))
      ii = ii+nExt(iSym)
    end if
  end do
  write(6,*)

  call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
  ip_Y = ip_X+nVV
  call FZero(Work(ip_X),nVV+nOA)
  call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,ip_X,ip_Y,.true.)
  call ChoMP2_Drv(irc,Dummy,CMO,EOcc,EVir)
  call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,ip_X,ip_Y,.false.)
  Wref = 0.0d0
  call ChoMP2_Drv(irc,EMP2,CMO,EOcc,EVir)
  if (irc /= 0) then
    write(6,*) 'LovMP2 failed'
    call Abend()
  end if
  iV = ip_X
  do iSym=1,nSym
    TrA(iSym) = ddot_(nExt(iSym),Work(iV),1+nExt(iSym),[1.0d0],0)
    iV = iV+nExt(iSym)**2
  end do
  call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
end if
write(6,*) '------------------------------------------------------'
write(6,*) ' Symm.   Tr(D):  Active        Frozen         Full    '
write(6,*) '------------------------------------------------------'
STrA = 0.0d0
STrF = 0.0d0
STrX = 0.0d0
do iSym=1,nSym
  write(6,'(2X,I4,9X,G11.4,3X,G11.4,3X,G11.4)') iSym,TrA(iSym),TrF(iSym),TrX(iSym)
  STrA = STrA+TrA(iSym)
  STrF = STrF+TrF(iSym)
  STrX = STrX+TrX(iSym)
end do
write(6,*) '------------------------------------------------------'
write(6,'(A,G11.4,3X,G11.4,3X,G11.4)') '          Sum: ',STrA,STrF,STrX
write(6,*) '------------------------------------------------------'
write(6,*)

write(6,'(A,F20.10,A)') ' Active region E2 contrib. = ',EMP2,' a.u.'
write(6,'(A,F20.10,A)') ' (Opposite-Spin contrib.   = ',-EOSMP2,' )'
write(6,*)

XEMP2 = EFRO
EMP2 = EMP2+EFRO
EOSMP2 = EOSMP2+EOSF

! Update runfile for subsequent calcs (e.g., CHCC)
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = 1+ioff
  ito = kEOcc+koff+nFro(iSym)
  call dcopy_(nOcc(iSym),EOcc(ifr),1,Work(ito),1)
  ifr = 1+joff
  ito = ito+nOcc(iSym)
  call dcopy_(nExt(iSym),EVir(ifr),1,Work(ito),1)
  ioff = ioff+nOcc(iSym)
  joff = joff+nExt(iSym)
  koff = koff+nBas(iSym)
end do
call Put_dArray('OrbE',Work(kEOcc),nxOrb)
call Put_dArray('Last orbitals',CMO,nSQ)

call GetMem('LCMO','Free','Real',iCMO,3*nSQ)
call GetMem('Saa','Free','Real',ipSaa,nxOrb)
call GetMem('Eorb','Free','Real',kEOcc,2*nxOrb)

2000 continue
if (min(iDo,jDo) == 0) then
  write(6,*)
  write(6,*) ' None of the occupied or virtual orbitals has been '
  write(6,*) ' assigned to the Active region of the molecule.    '
  write(6,*) ' This is presumably NOT what you want !!!          '
  write(6,*) ' MP2 will Stop here. Bye Bye !! '
  write(6,*)
  call Abend()
end if

call GetMem('SMAT','FREE','REAL',ipSQ,nSQ)

return

end subroutine LovMP2_Drv
