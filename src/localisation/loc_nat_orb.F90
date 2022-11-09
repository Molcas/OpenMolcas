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
!              by the user with the keyword LOCN.                      *
!              The threshold used for the orbitals splitting criterion *
!              is also required within the keyword LOCN.               *
!                                                                      *
!     Author: F. Aquilante   (Geneva, Feb. 2008)                       *
!                                                                      *
!***********************************************************************

use Localisation_globals, only: nActa, NamAct, BName, nBas, nFro, nSym, ThrSel
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: Cmo(*)
real(kind=wp), intent(inout) :: Xmo(*), OccN(*)
integer(kind=iwp), intent(in) :: mOrb(*)
integer(kind=iwp) :: i, ia, iab, iComp, ifr, iOff, iOpt, iSym, isymlbl, j, ja, jb, jC, jfr, jOcc, jOff, jto, jX, jZ, k, ka, kl, &
                     km, kOff, lScr, mOx, n_KO, n_OK, nBa, nBax, nBmx, nBx, nnB, nOrbmx, nOx
character(len=len(NamAct)) :: tmp
character(len=8) :: Label
integer(kind=iwp), allocatable :: jD(:), kD(:), lD(:)
real(kind=wp), allocatable :: C(:), CC(:), FOcc(:), S(:), Scr(:), SQ(:), Q(:), U(:), X(:), Z(:)
real(kind=wp), external :: ddot_

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
call mma_allocate(jD,nBmx,label='jD')
call mma_allocate(kD,nOrbmx,label='kD')
call mma_allocate(lD,nOrbmx,label='lD')
call mma_allocate(S,nnB,label='S')
call mma_allocate(SQ,nBmx**2,label='SQ')
isymlbl = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Label,iComp,S,isymlbl)
if (irc /= 0) return
lScr = nBmx*nOrbmx
call mma_allocate(C,lScr,label='C')
call mma_allocate(CC,lScr,label='CC')
call mma_allocate(X,lScr,label='X')
call mma_allocate(Z,lScr,label='Z')
call mma_allocate(Scr,lScr,label='Scr')
call mma_allocate(U,nOrbmx**2,label='U')
call mma_allocate(Q,nOrbmx,label='Q')
call mma_allocate(FOcc,nOrbmx,label='FOcc')

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
        nBa = nBa+1
        jD(nBa) = ia
      end if
    end do
  end do
  do ia=1,nBa
    ifr = jOff+jD(ia)
    call dcopy_(mOrb(iSym),Xmo(ifr),nBas(iSym),C(ia),nBa)
  end do
  do ia=1,nBa
    jb = jD(ia)
    jfr = kOff+jb*(jb-1)/2+1
    jto = nBas(iSym)*(ia-1)+1
    call dcopy_(jb,S(jfr),1,SQ(jto),1)
    jto = jto+jb
    do ka=jb+1,nBas(iSym)
      iab = kOff+ka*(ka-1)/2+jb
      SQ(jto) = S(iab)
      jto = jto+1
    end do
  end do
  nBx = max(1,nBas(iSym))
  nBax = max(1,nBa)
  call DGEMM_('T','N',nBa,mOrb(iSym),nBas(iSym),One,SQ,nBx,Xmo(jOff+1),nBx,Zero,Z,nBax)
  do i=1,mOrb(iSym)
    jC = nBa*(i-1)+1
    Q(i) = ddot_(nBa,C(jC),1,Z(jC),1)**2
  end do
  n_OK = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    jfr = jOff+nBas(iSym)*(i-1)+1
    if (sqrt(Q(i)) >= ThrSel) then
      jX = nBas(iSym)*n_OK+1
      call dcopy_(nBas(iSym),Xmo(jfr),1,X(jX),1)
      n_OK = n_OK+1
      kD(n_OK) = i
    else
      jZ = nBas(iSym)*n_KO+1
      call dcopy_(nBas(iSym),Xmo(jfr),1,Z(jZ),1)
      n_KO = n_KO+1
      lD(n_KO) = i
    end if
  end do
  call Square(S(kOff+1),SQ,1,nBas(iSym),nBas(iSym))
  mOx = max(1,mOrb(iSym))
  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,Cmo(jOff+1),nBx,SQ,nBx,Zero,CC,mOx)

  call DGEMM_('N','N',mOrb(iSym),n_OK,nBas(iSym),One,CC,mOx,X,nBx,Zero,U,mOx)
  jOcc = iOff+nFro(iSym)+1
  call Get_Nat_Lorb(OccN(jOcc),FOcc,n_OK,mOrb(iSym),kD,U)
  nOx = max(1,n_OK)
  call DGEMM_('N','N',nBas(iSym),n_OK,n_OK,One,X,nBx,U,nOx,Zero,Scr,nBx)
  do i=1,n_OK
    kl = nBas(iSym)*(i-1)+1
    km = jOff+nBas(iSym)*(kD(i)-1)+1
    call dcopy_(nBas(iSym),Scr(kl),1,Xmo(km),1)
  end do

  call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,CC,mOx,Z,nBx,Zero,U,mOx)
  call Get_Nat_Lorb(OccN(jOcc),FOcc,n_KO,mOrb(iSym),lD,U)
  nOx = max(1,n_KO)
  call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Z,nBx,U,nOx,Zero,Scr,nBx)
  do i=1,n_KO
    kl = nBas(iSym)*(i-1)+1
    j = lD(i)
    km = jOff+nBas(iSym)*(j-1)+1
    call dcopy_(nBas(iSym),Scr(kl),1,Xmo(km),1)
    k = jOcc+j-1
    OccN(k) = FOcc(j)
  end do
  do i=1,n_OK
    j = kD(i)
    k = jOcc+j-1
    OccN(k) = FOcc(j)
  end do

  iOff = iOff+nBas(iSym)
  jOff = jOff+nBas(iSym)*mOrb(iSym)
  kOff = kOff+nBas(iSym)*(nBas(iSym)+1)/2
end do

call mma_deallocate(jD)
call mma_deallocate(kD)
call mma_deallocate(lD)
call mma_deallocate(S)
call mma_deallocate(SQ)
call mma_deallocate(C)
call mma_deallocate(CC)
call mma_deallocate(X)
call mma_deallocate(Z)
call mma_deallocate(Scr)
call mma_deallocate(U)
call mma_deallocate(Q)
call mma_deallocate(FOcc)

return

end subroutine Loc_Nat_orb
