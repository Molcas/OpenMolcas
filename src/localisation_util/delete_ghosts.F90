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

subroutine Delete_Ghosts(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BName,nUniqAt,ThrS,isCASPT2,CMO,EOrb)
!***********************************************************************
!                                                                      *
! Purpose:  Eliminates MOs of ghost atoms from PT2 treatment           *
!                                                                      *
! Author:   F. Aquilante  (Geneva, July 2010)                          *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: irc, nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nDel(nSym), nUniqAt
character(len=LenIn8) :: BName(*)
real(kind=wp) :: ThrS, CMO(*), EOrb(*)
logical(kind=iwp) :: isCASPT2
#include "itmax.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ia, iAt, iBat, iC, iCMO, iD, ifr, ik, iOff, ip_iD, ip_nBas_per_Atom, ip_nBas_Start, ipAsh, ipQ, ipQa, ipS, &
                     ipSQ, ipZ, iQ, iQa, iQQ, iS, iSQ, iSym, isymlbl, ito, iv, iX, iZ, j, ja, jAt, jb, jBas, jBat, jC, jCMO, jfr, &
                     jjCMO, jjZ, jOff, jQ, jto, jX, jZ, kBas, kCMO, kOff, l_nBas_per_Atom, l_nBas_Start, lBas, LCMO, lsq, ltri, &
                     mAsh, n_KO, n_OK(8), nActa, nAk, nBa, nBasT, nBat, nBax, nBk, nBmx, nBx, NCMO, nOkk, nSmall, nSmx, nSQ, ntri
character(len=LenIn) :: NamAct(mxAtom), tmp
real(kind=r8), external :: ddot_

irc = 0

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nBasT = 0
ntri = 0
nSQ = 0
nBmx = 0
nSmx = 0
mAsh = 0
do i=1,nSym
  nBasT = nBasT+nBas(i)
  ntri = ntri+nBas(i)*(nBas(i)+1)/2
  nSQ = nSQ+nBas(i)**2
  nBmx = max(nBmx,nBas(i))
  nSmx = max(nSmx,nSsh(i))
  mAsh = max(mAsh,nIsh(i)+nAsh(i))
end do
NCMO = nSQ
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

! nUniqAt = # of symm. unique atoms. Initialize NamAct to blanks.
! ---------------------------------------------------------------

if ((nUniqAt < 1) .or. (nUniqAt > MxAtom)) then
  write(u6,'(A,I9)') 'nUniqAt =',nUniqAt
  call Abend()
end if
do iAt=1,nUniqAt
  NamAct(iAt) = ' '
end do

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

l_nBas_per_Atom = nUniqAt
l_nBas_Start = nUniqAt
call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,l_nBas_Start)

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
call GetMem('SMAT','ALLO','REAL',ipSQ,nSQ)
call GetMem('SLT','ALLO','REAL',ipS,nTri)
isymlbl = 1
call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
if (irc /= 0) return
ltri = 0
lsq = 0
do iSym=1,nSym
  call Square(Work(ipS+ltri),Work(ipSQ+lsq),1,nBas(iSym),nBas(iSym))
  ltri = ltri+nBas(iSym)*(nBas(iSym)+1)/2
  lsq = lsq+nBas(iSym)**2
end do
call GetMem('SLT','FREE','REAL',ipS,nTri)

call GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
call dcopy_(NCMO,CMO,1,WORK(LCMO),1)

!----------------------------------------------------------------------*
!     Compute Mulliken atomic charges of each occupied orbital         *
!             on each center to define the Active Site                 *
!----------------------------------------------------------------------*
call GetMem('Qai','Allo','Real',ipQ,nUniqAt*(mAsh+1))
ipQa = ipQ+nUniqAt*mAsh
call Fzero(Work(ipQa),nUniqAt)
call GetMem('Zm','Allo','Real',ipZ,nBmx*mAsh)
lBas = 0
iOff = 0
do iSym=1,nSym
  nOkk = nIsh(iSym)+nAsh(iSym)
  iSQ = ipSQ+iOff
  ipAsh = LCMO+iOff+nBas(iSym)*nFro(iSym)
  nBx = max(1,nBas(iSym))
  call DGEMM_('N','N',nBas(iSym),nOkk,nBas(iSym),One,Work(iSQ),nBx,Work(ipAsh),nBx,Zero,Work(ipZ),nBx)
  jBas = lBas+1
  kBas = lBas+nBas(iSym)
  call BasFun_Atom_(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),BName,jBas,kBas,nUniqAt,.false.)
  do ik=0,nOkk-1
    nAk = nUniqAt*ik
    nBk = nBas(iSym)*ik
    jCMO = ipAsh+nBk-1
    jZ = ipZ+nBk-1
    do iAt=0,nUniqAt-1
      iBat = iWork(ip_nBas_Start+iAt)
      jjCMO = jCMO+iBat
      jjZ = jZ+iBat
      nBat = iWork(ip_nBas_per_Atom+iAt)
      iQQ = ipQ+nAk+iAt
      Work(iQQ) = ddot_(nBat,Work(jjCMO),1,Work(jjZ),1)
    end do
  end do
  do iAt=0,nUniqAt-1
    jQ = ipQ+iAt
    iQa = ipQa+iAt
    Work(iQa) = Work(iQa)+ddot_(nOkk,Work(jQ),nUniqAt,Work(jQ),nUniqAt)
    if (sqrt(Work(iQa)) >= ThrS) then
      jBat = iWork(ip_nBas_Start+iAt)+lBas
      if (iWork(ip_nBas_per_Atom+iAt) > 0) NamAct(iAt+1) = BName(jBat)(1:LenIn)
    end if
  end do
  lBas = lBas+nBas(iSym)
  iOff = iOff+nBas(iSym)**2
end do
call GetMem('Zm','Free','Real',ipZ,nBmx*mAsh)
call GetMem('Qai','Free','Real',ipQ,nUniqAt*(mAsh+1))

!     We have now completed the definition of the active site
!----------------------------------------------------------------------*
call GetMem('ID_A','Allo','Inte',iD,nUniqAt)
nActa = 0
do iAt=1,nUniqAt
  if (NamAct(iAt)(1:4) /= ' ') then
    iWork(iD+nActa) = iAt
    nActa = nActa+1
  end if
end do
do iAt=1,nActa
  jAt = iWork(iD+iAt-1)
  NamAct(iAt) = NamAct(jAt)
end do
do iAt=nActa+1,nUniqAt
  NamAct(iAt)(1:4) = ' '
end do
write(u6,*)
write(u6,'(A,F6.3)') ' Threshold for atom selection: ',ThrS
write(u6,*)
if (nActa /= 0) then
  write(u6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
  write(u6,*)
  write(u6,*) (NamAct(i),i=1,nActa)
  write(u6,*)
else
  write(u6,*) ' None of the occupied non-frozen orbitals has been '
  write(u6,*) ' assigned to the Active region of the molecule.    '
  write(u6,*) ' This is presumably NOT what you want !!!          '
  write(u6,*) ' I will Stop here. Bye Bye !! '
  write(u6,*)
  call Abend()
end if

call GetMem('ID_A','Free','Inte',iD,nUniqAt)
call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)

!----------------------------------------------------------------------*
!     Virtual orbital selection                                        *
!----------------------------------------------------------------------*
call GetMem('ID_vir','Allo','Inte',ip_iD,nBmx+2*nSmx)
nSmall = nBmx**2+nSmx+3*nBmx*nSmx
call GetMem('Small','Allo','Real',iS,nSmall)
iQ = iS+nBmx**2
iC = iQ+nSmx
iZ = iC+nBmx*nSmx
iX = iZ+nBmx*nSmx
!
iOff = 0
kOff = 0
do iSym=1,nSym

  nBa = 0
  do ia=1,nBas(iSym)
    ja = ia+iOff
    tmp = BName(ja)(1:LenIn)
    do j=1,nActa
      if (NamAct(j) == tmp) then
        iWork(ip_iD+nBa) = ia
        nBa = nBa+1
      end if
    end do
  end do

  iCMO = LCMO+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
  do ia=1,nBa
    ifr = iCMO+iWork(ip_iD-1+ia)-1
    ito = iC+ia-1
    call dcopy_(nSsh(iSym),Work(ifr),nBas(iSym),Work(ito),nBa)
  end do

  iSQ = ipSQ+kOff
  do ia=1,nBa
    jb = iWork(ip_iD-1+ia)
    jfr = iSQ+nBas(iSym)*(jb-1)
    jto = iS+nBas(iSym)*(ia-1)
    call dcopy_(nBas(iSym),Work(jfr),1,Work(jto),1)
  end do

  nBx = max(1,nBas(iSym))
  nBax = max(1,nBa)
  call DGEMM_('T','N',nBa,nSsh(iSym),nBas(iSym),One,Work(iS),nBx,Work(iCMO),nBx,Zero,Work(iZ),nBax)
  do i=0,nSsh(iSym)-1
    jQ = iQ+i
    jC = iC+nBa*i
    jZ = iZ+nBa*i
    Work(jQ) = ddot_(nBa,Work(jC),1,Work(jZ),1)**2
  end do
  n_OK(iSym) = 0
  n_KO = 0
  do i=1,nSsh(iSym)
    jQ = iQ+i-1
    jfr = iCMO+nBas(iSym)*(i-1)
    if (sqrt(Work(jQ)) >= ThrS) then
      jX = iX+nBas(iSym)*n_OK(iSym)
      call dcopy_(nBas(iSym),Work(jfr),1,Work(jX),1)
      iWork(ip_iD+nBmx+n_OK(iSym)) = i
      n_OK(iSym) = n_OK(iSym)+1
    else
      jZ = iZ+nBas(iSym)*n_KO
      call dcopy_(nBas(iSym),Work(jfr),1,Work(jZ),1)
      iWork(ip_iD+nBmx+nSmx+n_KO) = i
      n_KO = n_KO+1
    end if
  end do

  call dcopy_(nBas(iSym)*n_OK(iSym),Work(iX),1,Work(iCMO),1)
  kCMO = iCMO+nBas(iSym)*n_OK(iSym)
  call dcopy_(nBas(iSym)*n_KO,Work(iZ),1,Work(kCMO),1)
  if (.not. isCASPT2) then
    jZ = iZ
    jOff = iOff+nFro(iSym)+nOkk
    do i=nBmx+1,nBmx+n_OK(iSym)
      iv = iWork(ip_iD-1+i)+jOff
      Work(jZ) = EOrb(iv)
      jZ = jZ+1
    end do
    do i=nBmx+nSmx+1,nBmx+nSmx+n_KO
      iv = iWork(ip_iD-1+i)+jOff
      Work(jZ) = EOrb(iv)
    end do
    call dcopy_(nSsh(iSym),Work(iZ),1,EOrb(1+jOff),1)
  end if

  iOff = iOff+nBas(iSym)
  kOff = kOff+nBas(iSym)**2
end do
call GetMem('Small','Free','Real',iS,nSmall)
call GetMem('ID_vir','Free','Inte',ip_iD,nBmx+2*nSmx)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Update nSsh, nDel for the Active site PT2
do iSym=1,nSym
  nDel(iSym) = nDel(iSym)+nSsh(iSym)-n_OK(iSym)
  nSsh(iSym) = n_OK(iSym)
end do

call dcopy_(NCMO,WORK(LCMO),1,CMO,1)

call GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
call GetMem('SMAT','FREE','REAL',ipSQ,nSQ)

return

end subroutine Delete_Ghosts
