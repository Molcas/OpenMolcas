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

subroutine BitMap_Localisation_Atom(PreFix)

use Localisation_globals, only: AnaNrm, ipCMO, ipMOrig, BName, nAtoms, nBas, nFro, nOrb2Loc, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
character(len=2), intent(in) :: PreFix
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iTyp, kC0(2), kC1, kC2, kX1
logical(kind=iwp) :: Debug
character(len=4) :: Typ(2)
character(len=12) :: BasNam
integer(kind=iwp), allocatable :: nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: Coord(:,:), CAt(:), DAt(:), Den(:), XAt(:)
character(len=24), parameter :: SecNam = 'BitMap_Localisation_Atom'

Debug = .false.

! Symmetry is not allowed!
! ------------------------

if (nSym /= 1) then
  call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
end if

! Allocate max. sym. block of density matrix
! and atom based density and CMO matrices.
! ------------------------------------------

call mma_allocate(Den,nBas(1)**2,label='BMpLoc')
call mma_allocate(DAt,nAtoms**2,label='DAt')
call mma_allocate(CAt,nAtoms*nOrb2Loc(1),label='CAt')
call mma_allocate(XAt,nAtoms*nOrb2Loc(1),label='XAt')

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtoms,label='nB_per_Atom')
call mma_allocate(nBas_Start,nAtoms,label='nB_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,BName,nBas(1),nAtoms,Debug)

! Compute density matrix, Den = CC^T, and set atom based matrices.
! Generate bitmap and perform sparsity analysis.
! ----------------------------------------------------------------

kC1 = ipMOrig+nBas(1)*nFro(1)
call GetDens_Localisation(Den,Work(kC1),nBas(1),nOrb2Loc(1))
call GetAt_Localisation(Den,nBas(1),nBas(1),DAt,nAtoms,2,nBas_per_Atom,nBas_Start,AnaNrm)
call GetAt_Localisation(Work(kC1),nBas(1),nOrb2Loc(1),CAt,nAtoms,1,nBas_per_Atom,nBas_Start,AnaNrm)
kX1 = ipCMO+nBas(1)*nFro(1)
call GetAt_Localisation(Work(kX1),nBas(1),nOrb2Loc(1),XAt,nAtoms,1,nBas_per_Atom,nBas_Start,AnaNrm)
call GenBMp_Localisation(DAt,CAt,XAt,nAtoms,1,'r','r','r',PreFix)
call Anasize_Localisation(DAt,CAt,XAt,nAtoms,nOrb2Loc(1),1)
write(u6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! Allocate memory for nuclear coordinates.
! Read nuclear coordinates from the runfile.
! Generate gnuplot files.
! ------------------------------------------

call mma_allocate(Coord,3,nAtoms,label='NucCoord')
call Get_dArray('Unique Coordinates',Coord,3*nAtoms)

call GetAt_Localisation(Den,nBas(1),nBas(1),DAt,nAtoms,2,nBas_per_Atom,nBas_Start,AnaNrm)
write(BasNam,'(A2,A10)') PreFix,'TotDensity'
call GenGnu_Localisation(BasNam,DAt,Coord,nAtoms)
Typ(1) = 'Dini'
Typ(2) = 'Dloc'
kC0(1) = ipMOrig
kC0(2) = ipCMO
do iTyp=1,2
  kC1 = kC0(iTyp)+nBas(1)*nFro(1)
  do i=1,nOrb2Loc(1)
    kC2 = kC1+nBas(1)*(i-1)
    call GetDens_Localisation(Den,Work(kC2),nBas(1),1)
    call GetAt_Localisation(Den,nBas(1),nBas(1),DAt,nAtoms,2,nBas_per_Atom,nBas_Start,AnaNrm)
    write(BasNam,'(A2,A4,I6)') PreFix,Typ(iTyp),i
    call GenGnu_Localisation(BasNam,DAt,Coord,nAtoms)
  end do
end do
write(u6,*) 'Gnuplot files have been generated. Norm: ',AnaNrm

! De-allocations.
! ---------------

call mma_deallocate(Den)
call mma_deallocate(CAt)
call mma_deallocate(DAt)
call mma_deallocate(XAt)
call mma_deallocate(nBas_per_Atom)
call mma_deallocate(nBas_Start)
call mma_deallocate(Coord)

end subroutine BitMap_Localisation_Atom
