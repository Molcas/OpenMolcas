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

subroutine Get_Orb_Select(irc,CMO,XMO,Eorb,Smat,Saa,BName,NamAct,nSym,nActa,mOrb,nBas,ortho,ThrSel,n_OK)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, r8

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nSym, nActa, mOrb(nSym), nBas(nSym)
integer(kind=iwp), intent(out) :: irc, n_OK(nSym)
real(kind=wp), intent(inout) :: CMO(*), Eorb(*)
real(kind=wp), intent(in) :: XMO(*), Smat(*), Saa(*), ThrSel
character(len=LenIn8), intent(in) :: BName(*)
character(len=LenIn), intent(in) :: NamAct(nActa)
logical(kind=iwp), intent(in) :: ortho
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ia, ifr, iOff, ip_C, ip_CC, ip_Fock, ip_iD, ip_U, ip_X, ipScr, iQ, iS, iSQ, iSym, ito, iZ, j, ja, jb, jC, &
                     jfr, jOff, jp_Fock, jQ, jto, jX, jZ, k, km, kOff, l, lOff, lScr, mOx, n_KO, nBa, nBax, nBmx, nBx, nORbmx, nOx
real(kind=wp) :: ThrS
character(len=LenIn) :: tmp
real(kind=r8), external :: ddot_

irc = 0

nBmx = 0
nOrbmx = 0
do iSym=1,nSym
  nBmx = max(nBmx,nBas(iSym))
  nOrbmx = max(nOrbmx,mOrb(iSym))
end do
call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx+nBmx)
call GetMem('Smx','Allo','Real',iSQ,nBmx**2)
call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+3)*nOrbmx)
lScr = nBmx*nOrbmx
ip_CC = ip_C+lScr
ip_X = ip_CC+lScr
iZ = ip_X+lScr
ipScr = iZ+lScr
ip_U = ipScr+lScr
iQ = ip_U+nOrbmx**2
ip_Fock = iQ+nOrbmx
jp_Fock = ip_Fock+nOrbmx
call Fzero(Work(iQ),nOrbmx)

iOff = 0
jOff = 0
kOff = 0
lOff = 0
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
  do ia=1,nBa
    ifr = jOff+iWork(ip_iD-1+ia)
    ito = ip_C+ia-1
    call dcopy_(mOrb(iSym),Xmo(ifr),nBas(iSym),Work(ito),nBa)
  end do
  iS = kOff+1
  do ia=1,nBa
    jb = iWork(ip_iD-1+ia)
    jfr = iS+nBas(iSym)*(jb-1)
    jto = iSQ+nBas(iSym)*(ia-1)
    call dcopy_(nbas(iSym),Smat(jfr),1,Work(jto),1)
  end do
  nBx = max(1,nBas(iSym))
  nBax = max(1,nBa)
  call DGEMM_('T','N',nBa,mOrb(iSym),nBas(iSym),One,Work(iSQ),nBx,Xmo(jOff+1),nBx,Zero,Work(iZ),nBax)
  do i=0,mOrb(iSym)-1
    jQ = iQ+i
    jC = ip_C+nBa*i
    jZ = iZ+nBa*i
    Work(jQ) = ddot_(nBa,Work(jC),1,Work(jZ),1)**2
  end do
  n_OK(iSym) = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    ThrS = ThrSel*Saa(lOff+i)
    jQ = iQ+i-1
    jfr = jOff+nBas(iSym)*(i-1)+1
    if (sqrt(Work(jQ)) >= ThrS) then
      jX = ip_X+nBas(iSym)*n_OK(iSym)
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
      iWork(ip_iD+nBmx+n_OK(iSym)) = i
      n_OK(iSym) = n_OK(iSym)+1
    else
      jZ = iZ+nBas(iSym)*n_KO
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
      iWork(ip_iD+nBmx+nOrbmx+n_KO) = i
      n_KO = n_KO+1
    end if
  end do

  mOx = max(1,mOrb(iSym))

  if (.not. ortho) then

    call Ortho_orb(Work(ip_X),Smat(iS),nBas(iSym),n_OK(iSym),2,.false.)
    call Ortho_orb(Work(iZ),Smat(iS),nBas(iSym),n_KO,2,.false.)
  end if

  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,Cmo(jOff+1),nBx,Smat(iS),nBx,Zero,Work(ip_CC),mOx)

  if (n_KO > 0) then
    call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,Work(ip_CC),mOx,Work(iZ),nBx,Zero,Work(ip_U),mOx)

    call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_KO,mOrb(iSym),iWork(ip_iD+nBmx+nOrbmx),Work(ip_U),iSym)

    call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Work(iZ),nBx,Work(ip_U),n_KO,Zero,Work(ipScr),nBx)

    ! Reorder the final MOs such that those of the active site come first
    km = jOff+nBas(iSym)*n_OK(iSym)+1
    call dcopy_(nBas(iSym)*n_KO,Work(ipScr),1,Cmo(km),1)
    call dcopy_(nOrbmx,Work(ip_Fock),1,Work(jp_Fock),1)
  end if

  call DGEMM_('N','N',mOrb(iSym),n_OK(iSym),nBas(iSym),One,Work(ip_CC),mOx,Work(ip_X),nBx,Zero,Work(ip_U),mOx)

  call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_OK(iSym),mOrb(iSym),iWork(ip_iD+nBmx),Work(ip_U),iSym)

  nOx = max(1,n_OK(iSym))
  call DGEMM_('N','N',nBas(iSym),n_OK(iSym),n_OK(iSym),One,Work(ip_X),nBx,Work(ip_U),nOx,Zero,Work(ipScr),nBx)

  do i=1,n_OK(iSym)
    j = iWork(ip_iD+nBmx-1+i)
    k = lOff+i
    l = ip_Fock+j-1
    Eorb(k) = Work(l)
  end do
  km = jOff+1
  call dcopy_(nBas(iSym)*n_OK(iSym),Work(ipScr),1,Cmo(km),1)

  do i=1,n_KO
    j = iWork(ip_iD+nBmx+nOrbmx-1+i)
    k = lOff+n_OK(iSym)+i
    l = jp_Fock+j-1
    Eorb(k) = Work(l)
  end do

  iOff = iOff+nBas(iSym)
  jOff = jOff+nBas(iSym)*mOrb(iSym)
  kOff = kOff+nBas(iSym)**2
  lOff = lOff+mOrb(iSym)
end do

call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+3)*nOrbmx)
call GetMem('Smx','Free','Real',iSQ,nBmx**2)
call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx+nBmx)

return

end subroutine Get_Orb_Select
