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

implicit real*8(a-h,o-z)
#include "files_mbpt2.fh"
#include "trafo.fh"
#include "corbinf.fh"
!defining One etc.
#include "real.fh"
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "chomp2_cfg.fh"

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

mAdOcc(1) = ipEOcc
do iSym=2,nSym
  mAdOcc(iSym) = mAdOcc(iSym-1)+nOcc(iSym-1)
end do

! The deleted "energies" are stached at the end!

mAdVir(1) = ipEVir
do iSym=2,nSym
  mAdVir(iSym) = mAdVir(iSym-1)+nExt(iSym-1)
end do

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
! Print occupied and virtual energies
do iSym=1,nSym
  if (nOcc(iSym) > 0) then
    write(6,*)
    write(6,*) 'Occupied energies, Sym',iSym
    write(6,*)
    do iOcc=0,nOcc(iSym)-1
      write(6,*) Work(mAdOcc(iSym)+iOcc)
    end do
  end if
  if (nExt(iSym) > 0) then
    write(6,*)
    write(6,*) 'Virtual energies, Sym',iSym
    write(6,*)
    do iExt=0,nExt(iSym)-1
      write(6,*) Work(mAdVir(iSym)+iExt)
    end do
  end if
  if (nFro(iSym) > 0) then
    write(6,*)
    write(6,*) 'Frozen energies, Sym',iSym
    write(6,*)
    do iFro=0,nFro(iSym)-1
      write(6,*) Work(mAdFro(iSym)+iFro)
    end do
  end if
  if (nDel(iSym) > 0) then
    write(6,*)
    write(6,*) 'Deleted energies, Sym',iSym
    write(6,*)
    do iDel=0,nDel(iSym)-1
      write(6,*) Work(mAdDel(iSym)+iDel)
      write(6,*)
    end do
  end if
end do
#endif

EMP2 = Zero
VECL2 = One

! Find the largest block needed to store exchange integrals
! either using IJ or AB as fix indeces.

nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do

lInt = nMaxOrb

! Allocate space for the block of exchange integrals
! type 1 and type 2 (ia|jb) and (ib|ja)

call GetMem('Int1','Allo','Real',ipInt1,lInt)
call GetMem('Int2','Allo','Real',ipInt2,lInt)

! Allocate a scratchvector for Coul and Exch subroutines.

call GetMem('Scr1','Allo','Real',ipScr1,lInt)

! Refer to the three-index storing of the symmetryblocks
! by this routine

#ifdef _DEBUGPRINT_
write(6,*) 'Before RHS_Mp2_help1'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),nI,nA)
  call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
end do
#endif
do iSymI=1,nSym
  do iSymJ=1,iSymI
    do iSymA=1,nSym
      iSymB = ieor(iSymA-1,(ieor(iSymI-1,iSymJ-1)))+1
      if (iSymB <= iSymA) then
        ! Check if there is any doubly excited configurations in
        ! the current symmetry.
        if ((nOrb(iSymI)+nDel(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))*(nOrb(iSymA)+nDel(iSymA))*(nOrb(iSymB)+nDel(iSymB)) /= 0) then
          call RHS_Mp2_help1(iSymA,iSymB,iSymI,iSymJ)
        end if !The constructed symmetry was empty
      end if   !No <occ occ | vir vir> in this symmetry
    end do     !ASym
  end do       !JSym
end do         !ISym
#ifdef _DEBUGPRINT_
write(6,*) 'After RHS_Mp2_help1'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),nI,nA)
  call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym
  do i=1,nOcc(iSym)+nExt(iSym)
    do j=1,nFro(iSym)
      Work(ip_Density(iSym)+i+nFro(iSym)-1+(j-1)*(nOrb(iSym)+nDel(iSym))) = Work(ip_Density(iSym)+j-1+(i+nFro(iSym)-1)* &
                                                                                 (nOrb(iSym)+nDel(iSym)))
    end do
    do j=1,nDel(iSym)
      Work(ip_Density(iSym)+j+nOrb(iSym)-1+(i+nFro(iSym)-1)*(nOrb(iSym)+nDel(iSym))) = Work(ip_Density(iSym)+i+nFro(iSym)-1+ &
                                                                                            (j+nOrb(iSym)-1)*(nOrb(iSym)+ &
                                                                                            nDel(iSym)))
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(6,*) 'After RHS_Mp2_help1 xx'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),nI,nA)
  call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
end do
#endif

! Since some terms in the Lagrangian is dependent of Pab and Pij we
! have to loop through all the symmetries again since they are not done
! until the loop is over.

do iSymI=1,nSym
  do iSymJ=1,iSymI
    do iSymA=1,nSym
      iSymB = ieor(iSymA-1,(ieor(iSymI-1,iSymJ-1)))+1
      if (iSymB <= iSymA) then
        ! Check if there is any doubly excited configurations in
        ! the current symmetry.
        if ((nOrb(iSymI)+nDel(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))*(nOrb(iSymA)+nDel(iSymA))*(nOrb(iSymB)+nDel(iSymB)) /= 0) then
          call RHS_Mp2_help2(iSymA,iSymB,iSymI,iSymJ)
        end if !The constructed symmetry was empty
      end if   !No <occ occ | vir vir> in this symmetry
    end do     !ASym
  end do       !JSym
end do         !ISym
#ifdef _DEBUGPRINT_
write(6,*) 'After RHS_Mp2_help1 yy'
do iSym=1,nSym
  nI = nOcc(iSym)+nFro(iSym)
  nA = nExt(iSym)+nDel(iSym)
  nB = nOrb(iSym)+nDel(iSym)
  call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),nI,nA)
  call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
end do
#endif

! Need two extra integral fields due to symmetrization

if (.not. NoGamma) then
  call GetMem('Int1_2','Allo','Real',ipInt1_2,lInt)
  call GetMem('Int2_2','Allo','Real',ipInt2_2,lInt)

  ! Construct and backtransform nonsep 2pdm.

  call Gamma_new()

  call GetMem('Int1_2','Free','Real',ipInt1_2,lInt)
  call GetMem('Int2_2','Free','Real',ipInt2_2,lInt)
end if

call GetMem('Int1','Free','Real',ipInt1,lInt)
call GetMem('Int2','Free','Real',ipInt2,lInt)
call GetMem('Scr1','Free','Real',ipScr1,lInt)

#ifdef _DEBUGPRINT_
write(6,*) 'EMP2 is ',EMP2
write(6,*) ' '
do iSym=1,nSym
  write(6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Work(ip_Density(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(6,*) 'WDensity matrix for Symm:',iSym
  call RecPrt('MP2WDensity','',Work(ip_WDensity(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(6,*) 'Lagrangian matrix for symm',iSym
  call RecPrt('Lagr2','',Work(ip_Mp2Lagr(iSym)),nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
end do
#endif

VECL2 = sqrt(One/VECL2)

return

end subroutine RHS_MP2
