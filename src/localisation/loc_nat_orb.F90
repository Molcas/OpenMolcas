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

subroutine Loc_Nat_orb(irc,Cmo,Xmo,OccN,mOrb)
!***********************************************************************
!                                                                      *
!     Purpose: compute Localized Natural Orbitals (LNOs) to be used    *
!              for instance in Effective Bond Order (EBO) analysis     *
!                                                                      *
!              The density matrix in localized MO basis is             *
!              diagonalized separately in 2 subblocks defined by the   *
!              orbitals that extend within two distinct subregions     *
!              (atoms) of the molecule.                                *
!              The gross Mulliken population of each orbital on the    *
!              atoms of the two subregions determines the splitting    *
!              of the orbitals.                                        *
!              The atoms defining the "active subregion" are specified *
!              by the user with the keyword LOCN .                     *
!              The threshold used for the orbitals splitting criterium *
!              is also required within the keyword LOCN .              *
!                                                                      *
!     Author: F. Aquilante   (Geneva, Feb. 2008)                       *
!                                                                      *
!***********************************************************************

use Localisation_globals, only: nActa, NamAct, BName, nBas, nFro, nSym, ThrSel
use Constants, only: Zero, One
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: Cmo(*)
real(kind=wp), intent(inout) :: Xmo(*), OccN(*)
integer(kind=iwp), intent(in) :: mOrb(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ia, iab, ifr, iOff, ip_C, ip_CC, ip_FOcc, ip_iD, ip_U, ip_X, ipS, ipScr, iQ, iS, iSQ, iSym, isymlbl, ito, &
                     iZ, j, ja, jb, jC, jfr, jOcc, jOff, jQ, jto, jX, jZ, k, ka, kl, km, kOff, l, lScr, mOx, n_KO, n_OK, nBa, &
                     nBax, nBmx, nBx, nnB, nOrbmx, nOx
character(len=len(NamAct)) :: tmp
real(kind=r8), external :: ddot_
!***********************************************************************
integer(kind=iwp) :: jD, kD, lD
jD(i) = iWork(ip_iD-1+i)
!*****
kD(i) = iWork(ip_iD+nBmx-1+i)
!*****
lD(i) = iWork(ip_iD+nBmx+nOrbmx-1+i)
!***********************************************************************

irc = 0

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
nnB = 0
nBmx = 0
nOrbmx = 0
do iSym=1,nSym
  nnB = nnB+nBas(iSym)*(nBas(iSym)+1)/2
  nBmx = max(nBmx,nBas(iSym))
  nOrbmx = max(nOrbmx,mOrb(iSym))
end do
call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx+nBmx)
call GetMem('Smat','Allo','Real',ipS,nnB+nBmx**2)
iSQ = ipS+nnB
isymlbl = 1
call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
if (irc /= 0) return
call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
lScr = nBmx*nOrbmx
ip_CC = ip_C+lScr
ip_X = ip_CC+lScr
iZ = ip_X+lScr
ipScr = iZ+lScr
ip_U = ipScr+lScr
iQ = ip_U+nOrbmx**2
ip_FOcc = iQ+nOrbmx
!
iOff = 0
jOff = 0
kOff = 0
do iSym=1,nSym

  nBa = 0
  do ia=1,nBas(iSym)
    ja = ia+iOff
    tmp = BName(ja)(1:len(tmp))
    do j=1,nActa
      if (NamAct(j) == tmp) then
        iWork(ip_iD+nBa) = ia
        nBa = nBa+1
      end if
    end do
  end do
  do ia=1,nBa
    ifr = jOff+jD(ia)
    ito = ip_C+ia-1
    call dcopy_(mOrb(iSym),Xmo(ifr),nBas(iSym),Work(ito),nBa)
  end do
  iS = ipS+kOff
  do ia=1,nBa
    jb = jD(ia)
    jfr = iS+jb*(jb-1)/2
    jto = iSQ+nBas(iSym)*(ia-1)
    call dcopy_(jb,Work(jfr),1,Work(jto),1)
    jto = jto+jb
    do ka=jb+1,nBas(iSym)
      iab = iS+ka*(ka-1)/2+jb-1
      Work(jto) = Work(iab)
      jto = jto+1
    end do
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
  n_OK = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    jQ = iQ+i-1
    jfr = jOff+nBas(iSym)*(i-1)+1
    if (sqrt(Work(jQ)) >= ThrSel) then
      jX = ip_X+nBas(iSym)*n_OK
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
      iWork(ip_iD+nBmx+n_OK) = i
      n_OK = n_OK+1
    else
      jZ = iZ+nBas(iSym)*n_KO
      call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
      iWork(ip_iD+nBmx+nOrbmx+n_KO) = i
      n_KO = n_KO+1
    end if
  end do
  call Square(Work(iS),Work(iSQ),1,nBas(iSym),nBas(iSym))
  mOx = max(1,mOrb(iSym))
  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,Cmo(jOff+1),nBx,Work(iSQ),nBx,Zero,Work(ip_CC),mOx)

  call DGEMM_('N','N',mOrb(iSym),n_OK,nBas(iSym),One,Work(ip_CC),mOx,Work(ip_X),nBx,Zero,Work(ip_U),mOx)
  jOcc = iOff+nFro(iSym)+1
  call Get_Nat_Lorb(OccN(jOcc),Work(ip_FOcc),n_OK,mOrb(iSym),iWork(ip_iD+nBmx),Work(ip_U),iSym)
  nOx = max(1,n_OK)
  call DGEMM_('N','N',nBas(iSym),n_OK,n_OK,One,Work(ip_X),nBx,Work(ip_U),nOx,Zero,Work(ipScr),nBx)
  do i=1,n_OK
    kl = ipScr+nBas(iSym)*(i-1)
    j = kD(i)
    km = jOff+nBas(iSym)*(j-1)+1
    call dcopy_(nBas(iSym),Work(kl),1,Xmo(km),1)
  end do

  call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,Work(ip_CC),mOx,Work(iZ),nBx,Zero,Work(ip_U),mOx)
  call Get_Nat_Lorb(OccN(jOcc),Work(ip_FOcc),n_KO,mOrb(iSym),iWork(ip_iD+nBmx+nOrbmx),Work(ip_U),iSym)
  nOx = max(1,n_KO)
  call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Work(iZ),nBx,Work(ip_U),nOx,Zero,Work(ipScr),nBx)
  do i=1,n_KO
    kl = ipScr+nBas(iSym)*(i-1)
    j = lD(i)
    km = jOff+nBas(iSym)*(j-1)+1
    call dcopy_(nBas(iSym),Work(kl),1,Xmo(km),1)
    k = jOcc+j-1
    l = ip_FOcc+j-1
    OccN(k) = Work(l)
  end do
  do i=1,n_OK
    j = kD(i)
    k = jOcc+j-1
    l = ip_FOcc+j-1
    OccN(k) = Work(l)
  end do

  iOff = iOff+nBas(iSym)
  jOff = jOff+nBas(iSym)*mOrb(iSym)
  kOff = kOff+nBas(iSym)*(nBas(iSym)+1)/2
end do

call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
call GetMem('Smat','Free','Real',ipS,nnB+nBmx**2)
call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx+nBmx)

return

end subroutine Loc_Nat_orb
