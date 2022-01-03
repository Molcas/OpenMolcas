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

subroutine Get_Vir_Select(irc,CMO,XMO,Eorb,Smat,BName,NamAct,ind_V,nSym,nActa,mOrb,nBas,ortho,n_OK)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit None
#include "Molcas.fh"
integer(kind=iwp) :: irc, ind_V(*), nSym, nActa, mOrb(nSym), nBas(nSym), n_OK(nSym)
real(kind=wp) :: CMO(*), XMO(*), Eorb(*), Smat(*)
character(len=LenIn8) :: BName(*)
character(len=LenIn) :: NamAct(nActa)
logical(kind=iwp) :: ortho
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iOff, ip_C, ip_CC, ip_Fock, ip_iD, ip_U, ip_X, ipScr, iS, iSym, iZ, j, ja, jfr, jOff, jp_Fock, jX, jZ, k, &
                     km, kOff, kx, l, lOff, lScr, mOx, n_KO, nBmx, nBx, nOrbmx, nOx
character(len=LenIn) :: tmp

irc = 0

nBmx = 0
nOrbmx = 0
do iSym=1,nSym
  nBmx = max(nBmx,nBas(iSym))
  nOrbmx = max(nOrbmx,mOrb(iSym))
end do
call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx)
call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
lScr = nBmx*nOrbmx
ip_CC = ip_C+lScr
ip_X = ip_CC+lScr
iZ = ip_X+lScr
ipScr = iZ+lScr
ip_U = ipScr+lScr
ip_Fock = ip_U+nOrbmx**2
jp_Fock = ip_Fock+nOrbmx

iOff = 0
jOff = 0
kOff = 0
lOff = 0
do iSym=1,nSym

  iS = kOff+1
  n_OK(iSym) = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    ja = iOff+ind_V(i+iOff)
    tmp = BName(ja)(1:LenIn)
    jfr = jOff+nBas(iSym)*(i-1)+1
    kx = 0

    !write(u6,*) ' We simulate Afreeze with all Vir'
    !kx = 1

    do j=1,nActa
      if (NamAct(j) == tmp) kx = kx+1
    end do
    if (kx > 0) then
      jX = ip_X+nBas(iSym)*n_OK(iSym)
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
      iWork(ip_iD+n_OK(iSym)) = i
      n_OK(iSym) = n_OK(iSym)+1
    else
      jZ = iZ+nBas(iSym)*n_KO
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
      iWork(ip_iD+nOrbmx+n_KO) = i
      n_KO = n_KO+1
    end if
  end do

  if (.not. ortho) then

    call Ortho_orb(Work(ip_X),Smat(iS),nBas(iSym),n_OK(iSym),2,.false.)
    call Ortho_orb(Work(iZ),Smat(iS),nBas(iSym),n_KO,2,.false.)
  end if

  nBx = max(1,nBas(iSym))
  mOx = max(1,mOrb(iSym))

  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,Cmo(jOff+1),nBx,Smat(iS),nBx,Zero,Work(ip_CC),mOx)

  if (n_KO > 0) then
    call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,Work(ip_CC),mOx,Work(iZ),nBx,Zero,Work(ip_U),mOx)

    call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_KO,mOrb(iSym),iWork(ip_iD+nOrbmx),Work(ip_U),iSym)

    call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Work(iZ),nBx,Work(ip_U),n_KO,Zero,Work(ipScr),nBx)

    ! Reorder the final MOs such that those of the active site come first
    km = jOff+nBas(iSym)*n_OK(iSym)+1
    call dcopy_(nBas(iSym)*n_KO,Work(ipScr),1,Cmo(km),1)
    call dcopy_(nOrbmx,Work(ip_Fock),1,Work(jp_Fock),1)
  end if

  call DGEMM_('N','N',mOrb(iSym),n_OK(iSym),nBas(iSym),One,Work(ip_CC),mOx,Work(ip_X),nBx,Zero,Work(ip_U),mOx)

  call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_OK(iSym),mOrb(iSym),iWork(ip_iD),Work(ip_U),iSym)

  nOx = max(1,n_OK(iSym))
  call DGEMM_('N','N',nBas(iSym),n_OK(iSym),n_OK(iSym),One,Work(ip_X),nBx,Work(ip_U),nOx,Zero,Work(ipScr),nBx)

  do i=1,n_OK(iSym)
    j = iWork(ip_iD-1+i)
    k = lOff+i
    l = ip_Fock+j-1
    Eorb(k) = Work(l)
  end do
  km = jOff+1
  call dcopy_(nBas(iSym)*n_OK(iSym),Work(ipScr),1,Cmo(km),1)

  do i=1,n_KO
    j = iWork(ip_iD+nOrbmx-1+i)
    k = lOff+n_OK(iSym)+i
    l = jp_Fock+j-1
    Eorb(k) = Work(l)
  end do

  iOff = iOff+nBas(iSym)
  jOff = jOff+nBas(iSym)*mOrb(iSym)
  kOff = kOff+nBas(iSym)**2
  lOff = lOff+mOrb(iSym)
end do

call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx)

return

end subroutine Get_Vir_Select
