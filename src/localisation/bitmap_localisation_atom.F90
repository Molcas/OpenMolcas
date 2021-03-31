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

implicit real*8(a-h,o-z)
character*2 PreFix
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

character*24 SecNam
parameter(SecNam='BitMap_Localisation_Atom')

character*4 Typ(2)
character*12 BasNam

integer kC0(2)

logical Debug

Debug = .false.

! Symmetry is not allowed!
! ------------------------

if (nSym /= 1) then
  call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
end if

! Allocate max. sym. block of density matrix
! and atom based density and CMO matrices.
! ------------------------------------------

lDen = nBas(1)**2
lDAt = nAtoms**2
lCAt = nAtoms*nOrb2Loc(1)
lXAt = lCAt
call GetMem('BMpLoc','Allo','Real',ipDen,lDen)
call GetMem('DAt','Allo','Real',ipDAt,lDAt)
call GetMem('CAt','Allo','Real',ipCAt,lCAt)
call GetMem('XAt','Allo','Real',ipXAt,lXAt)

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

l_nBas_per_Atom = nAtoms
l_nBas_Start = nAtoms
call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,l_nBas_Start)
call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),Name,nBas(1),nAtoms,Debug)

! Compute density matrix, Den = CC^T, and set atom based matrices.
! Generate bitmap and perform sparsity analysis.
! ----------------------------------------------------------------

kC1 = ipMOrig+nBas(1)*nFro(1)
call GetDens_Localisation(Work(ipDen),Work(kC1),nBas(1),nOrb2Loc(1))
call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),Work(ipDAt),nAtoms,2,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
call GetAt_Localisation(Work(kC1),nBas(1),nOrb2Loc(1),Work(ipCAt),nAtoms,1,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
kX1 = ipCMO+nBas(1)*nFro(1)
call GetAt_Localisation(Work(kX1),nBas(1),nOrb2Loc(1),Work(ipXAt),nAtoms,1,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
call GenBMp_Localisation(Work(ipDAt),Work(ipCAt),Work(ipXAt),nAtoms,1,'r','r','r',PreFix)
call Anasize_Localisation(Work(ipDAt),Work(ipCAt),Work(ipXAt),nAtoms,nOrb2Loc(1),1)
write(6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! Allocate memory for nuclear coordinates.
! Read nuclear coordinates from the runfile.
! Generate gnuplot files.
! ------------------------------------------

lCoord = 3*nAtoms
call GetMem('NucCoord','Allo','Real',ipCoord,lCoord)
call Get_dArray('Unique Coordinates',Work(ipCoord),lCoord)

call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),Work(ipDAt),nAtoms,2,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
write(BasNam,'(A2,A10)') PreFix,'TotDensity'
call GenGnu_Localisation(BasNam,Work(ipDAt),Work(ipCoord),nAtoms)
Typ(1) = 'Dini'
Typ(2) = 'Dloc'
kC0(1) = ipMOrig
kC0(2) = ipCMO
do iTyp=1,2
  kC1 = kC0(iTyp)+nBas(1)*nFro(1)
  do i=1,nOrb2Loc(1)
    kC2 = kC1+nBas(1)*(i-1)
    call GetDens_Localisation(Work(ipDen),Work(kC2),nBas(1),1)
    call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),Work(ipDAt),nAtoms,2,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
    write(BasNam,'(A2,A4,I6)') PreFix,Typ(iTyp),i
    call GenGnu_Localisation(BasNam,Work(ipDAt),Work(ipCoord),nAtoms)
  end do
end do
write(6,*) 'Gnuplot files have been generated. Norm: ',AnaNrm

! De-allocations.
! ---------------

call GetMem('NucCoord','Free','Real',ipCoord,lCoord)
call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)
call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('XAt','Free','Real',ipXAt,lXAt)
call GetMem('CAt','Free','Real',ipCAt,lCAt)
call GetMem('DAt','Free','Real',ipDAt,lDAt)
call GetMem('BMpLoc','Free','Real',ipDen,lDen)

end subroutine BitMap_Localisation_Atom
