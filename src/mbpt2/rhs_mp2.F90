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

subroutine RHS_MP2()
! The RHS for the MP2-gradients

use MBPT2_Global, only: Density, EMP2, LuIntM, mAdOcc, mAdVir, VECL2
use ChoMP2, only: NoGamma
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use MBPT2_Global, only: EOcc, EVir, mAdDel, mAdFro, Mp2Lagr, WDensity
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: i, iSym, iSym1, iSym2, iSymA, iSymB, iSymI, iSymJ, j, LIADOUT, nDelTot, nMaxOrb, nVirTot
real(kind=wp), allocatable :: Int1(:), Int1_2(:), Int2(:), Int2_2(:), Scr1(:)
#include "trafo.fh"
#include "corbinf.fh"

IAD13 = 0
LIADOUT = 3*36*36

! Read transformed integrals from a logical unit to a buffer

call iDAFILE(LuIntM,2,IADOUT,LIADOUT,IAD13)

! Sort startadresses for Virtual and occupied orbital
! energies in nice vectors

nVirTot = 0
nDelTot = 0
do i=1,nSym
  nVirTot = nVirTot+nExt(i)
  nDelTot = nDelTot+nDel(i)
end do

! The frozen energies are stashed at the end!

mAdOcc(1) = 1
do iSym=2,nSym
  mAdOcc(iSym) = mAdOcc(iSym-1)+nOcc(iSym-1)
end do

! The deleted "energies" are stashed at the end!

mAdVir(1) = 1
do iSym=2,nSym
  mAdVir(iSym) = mAdVir(iSym-1)+nExt(iSym-1)
end do

#ifdef _DEBUGPRINT_
! Print occupied and virtual energies
do iSym=1,nSym
  if (nOcc(iSym) > 0) then
    write(u6,*)
    write(u6,*) 'Occupied energies, Sym',iSym
    write(u6,*)
    do iOcc=0,nOcc(iSym)-1
      write(u6,*) EOcc(mAdOcc(iSym)+iOcc)
    end do
  end if
  if (nExt(iSym) > 0) then
    write(u6,*)
    write(u6,*) 'Virtual energies, Sym',iSym
    write(u6,*)
    do iExt=0,nExt(iSym)-1
      write(u6,*) EVir(mAdVir(iSym)+iExt)
    end do
  end if
  if (nFro(iSym) > 0) then
    write(u6,*)
    write(u6,*) 'Frozen energies, Sym',iSym
    write(u6,*)
    do iFro=0,nFro(iSym)-1
      write(u6,*) EOcc(mAdFro(iSym)+iFro)
    end do
  end if
  if (nDel(iSym) > 0) then
    write(u6,*)
    write(u6,*) 'Deleted energies, Sym',iSym
    write(u6,*)
    do iDel=0,nDel(iSym)-1
      write(u6,*) EVir(mAdDel(iSym)+iDel)
      write(u6,*)
    end do
  end if
end do
#endif

EMP2 = Zero
VECL2 = One

! Find the largest block needed to store exchange integrals
! either using IJ or AB as fix indices.

nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do

! Allocate space for the block of exchange integrals
! type 1 and type 2 (ia|jb) and (ib|ja)

call mma_allocate(Int1,nMaxOrb,label='Int1')
call mma_allocate(Int2,nMaxOrb,label='Int2')

! Allocate a scratchvector for Coul and Exch subroutines.

call mma_allocate(Scr1,nMaxOrb,label='Scr1')

! Refer to the three-index storing of the symmetryblocks
! by this routine

#ifdef _DEBUGPRINT_
write(u6,*) 'Before RHS_Mp2_help1'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Mp2Lagr%SB(iSym)%A1,nI,nA)
  call RecPrt('Dens',' ',Density%SB(iSym)%A1,nB,nB)
end do
#endif
do iSymI=1,nSym
  do iSymJ=1,iSymI
    do iSymA=1,nSym
      iSymB = Mul(iSymA,Mul(iSymI,iSymJ))
      if (iSymB <= iSymA) then
        ! Check if there is any doubly excited configurations in
        ! the current symmetry.
        if ((nOrb(iSymI)+nDel(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))*(nOrb(iSymA)+nDel(iSymA))*(nOrb(iSymB)+nDel(iSymB)) /= 0) then
          call RHS_Mp2_help1(iSymA,iSymB,iSymI,iSymJ,Int1,Int2,Scr1)
        end if !The constructed symmetry was empty
      end if   !No <occ occ | vir vir> in this symmetry
    end do     !ASym
  end do       !JSym
end do         !ISym
#ifdef _DEBUGPRINT_
write(u6,*) 'After RHS_Mp2_help1'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Mp2Lagr%SB(iSym)%A1,nI,nA)
  call RecPrt('Dens',' ',Density%SB(iSym)%A1,nB,nB)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym
  do i=1,nOcc(iSym)+nExt(iSym)
    do j=1,nFro(iSym)
      Density%SB(iSym)%A2(i+nFro(iSym),j) = Density%SB(iSym)%A2(j,i+nFro(iSym))
    end do
    do j=1,nDel(iSym)
      Density%SB(iSym)%A2(j+nOrb(iSym),i+nFro(iSym)) = Density%SB(iSym)%A2(i+nFro(iSym),j+nOrb(iSym))
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'After RHS_Mp2_help1 xx'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Mp2Lagr%SB(iSym)%A1,nI,nA)
  call RecPrt('Dens',' ',Density%SB(iSym)%A1,nB,nB)
end do
#endif

! Since some terms in the Lagrangian is dependent of Pab and Pij we
! have to loop through all the symmetries again since they are not done
! until the loop is over.

do iSymI=1,nSym
  do iSymJ=1,iSymI
    do iSymA=1,nSym
      iSymB = Mul(iSymA,Mul(iSymI,iSymJ))
      if (iSymB <= iSymA) then
        ! Check if there is any doubly excited configurations in
        ! the current symmetry.
        if ((nOrb(iSymI)+nDel(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))*(nOrb(iSymA)+nDel(iSymA))*(nOrb(iSymB)+nDel(iSymB)) /= 0) then
          call RHS_Mp2_help2(iSymA,iSymB,iSymI,iSymJ,Int1,Int2,Scr1)
        end if !The constructed symmetry was empty
      end if   !No <occ occ | vir vir> in this symmetry
    end do     !ASym
  end do       !JSym
end do         !ISym
#ifdef _DEBUGPRINT_
write(u6,*) 'After RHS_Mp2_help1 yy'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Mp2Lagr%SB(iSym)%A1,nI,nA)
  call RecPrt('Dens',' ',Density%SB(iSym)%A1,nB,nB)
end do
#endif

! Need two extra integral fields due to symmetrization

if (.not. NoGamma) then
  call mma_allocate(Int1_2,nMaxOrb,label='Int1_2')
  call mma_allocate(Int2_2,nMaxOrb,label='Int2_2')

  ! Construct and backtransform nonsep 2pdm.

  call Gamma_new(Int1,Int2,Int1_2,Int2_2,Scr1)

  call mma_deallocate(Int1_2)
  call mma_deallocate(Int2_2)
end if

call mma_deallocate(Int1)
call mma_deallocate(Int2)
call mma_deallocate(Scr1)

#ifdef _DEBUGPRINT_
write(u6,*) 'EMP2 is ',EMP2
write(u6,*) ' '
do iSym=1,nSym
  write(u6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Density%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(u6,*) 'WDensity matrix for Symm:',iSym
  call RecPrt('MP2WDensity','',WDensity%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(u6,*) 'Lagrangian matrix for symm',iSym
  call RecPrt('Lagr2','',Mp2Lagr%SB(iSyM)%A1,nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
end do
#endif

VECL2 = sqrt(One/VECL2)

return

end subroutine RHS_MP2
