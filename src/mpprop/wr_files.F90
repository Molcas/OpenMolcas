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

! OBSERVE! If the output to the mpprop file that this subroutine generates
!   is modified, then notify the person responsible for qmstat, since that
!   program uses the mpprop outputfile.
subroutine Wr_Files(nAtoms,nCenters,nMltPl,NORBI,NOCOB,NOCOB_b,OENE,OENE_b,LAllCenters)
!EB subroutine Wr_Files(nAtoms,nCenters,nMltPl,NORBI,NOCOB,OENE,iBond,

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nCenters, nMltPl, NORBI, NOCOB, NOCOB_b !, iBond(2,nCenters)
real(kind=wp), intent(in) :: OENE(NOCOB), OENE_b(NOCOB_b)
logical(kind=iwp), intent(in) :: LAllCenters
#include "MpData.fh"
#include "WrkSpc.fh"
#include "MolProp.fh"
integer(kind=iwp) :: i, iComp, iMltPl, j, k, Lu, nComp
real(kind=wp) :: MolPol(6)
character(len=8) :: FileName
logical(kind=iwp) :: Exists

do i=1,6
  MolPol(i) = Zero
end do

Lu = 30
FileName = 'MPPROP'
call OpnFl(FileName,Lu,Exists)
rewind(Lu)

write(Lu,6) '**************************************************'
write(Lu,6) '* Molecule'
write(Lu,6) Title
write(Lu,6) '* Method'
write(Lu,6) Method
write(Lu,6) '* Level of Multipoles and Polarizabilities'
write(Lu,3) nMltPl,1

! WE NOW HAVE THE MULTIPOLES IN TWO FORMS
!
! A) SUMMED ONTO ATOMIC TERMS ONLY
! B) SUMMED ONTO ATOMS AND BONDS
!
! EXPANSION TO BE USED FOR ELECTROSTATICS
!
!       CENTER
!       MULTIPOLE
!       POLARIZABILITY

if (.not. LAllCenters) then
  write(Lu,'(A)') '* Atomic centers '
  write(Lu,2) nAtoms
  do i=1,nAtoms
    write(Lu,3) iAtomType(i),iAtomPar(i),Cen_Lab(i*(i+1)/2)
    write(Lu,1) cor(1,i,i),cor(2,i,i),cor(3,i,i)
    do iMltPl=0,nMltPl
      nComp = (iMltPl+1)*(iMltPl+2)/2
      write(Lu,1) (Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+i-1),iComp=1,nComp)
    end do
    do iMltPl=0,2
      nComp = (iMltPl+1)*(iMltPl+2)/2
      write(Lu,1) (Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*(iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
    end do
    write(Lu,1) (Work(iAtPolAd+nAtoms*j+i-1),j=0,5)
    do j=0,5
      MolPol(j+1) = MolPol(j+1)+Work(iAtPolAd+nAtoms*j+i-1)
    end do
  end do

else

  write(Lu,6) '* All centers'
  write(Lu,2) nCenters
  do i=1,nAtoms
    write(Lu,3) iAtomType(i),iAtomPar(i),Cen_Lab(i*(i+1)/2)
    write(Lu,1) Cor(1,i,i),Cor(2,i,i),Cor(3,i,i)
    do iMltPl=0,nMltPl
      nComp = (iMltPl+1)*(iMltPl+2)/2
      write(Lu,1) (Work(iAtBoMltPlAd(iMltPl)+nCenters*(iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
    end do
    do iMltPl=0,2
      nComp = (iMltPl+1)*(iMltPl+2)/2
      write(Lu,1) (Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*(iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
    end do
    ! Begin EB
    ! if "pointer" iAtBoPolAd is not associated with memory
    ! print zero
    if (iAtBoPolAd == 0) then
      write(Lu,1) Zero,Zero,Zero,Zero,Zero,Zero
    else
    ! End EB
      write(Lu,1) (Work(iAtBoPolAd+nCenters*k+i*(i-1)/2+i-1),k=0,5)
    ! Begin EB
    end if
    ! End EB
  end do
  do i=1,nAtoms
    do j=1,i-1
      if (BondMat(i,j)) then
        write(Lu,5) 0,iBondPar(i),iAtomType(i),iAtomPar(i),iAtomType(j),iAtomPar(j),Cen_Lab(i*(i-1)/2+j)
        write(Lu,1) Cor(1,i,j),Cor(2,i,j),Cor(3,i,j)
        do iMltPl=0,nMltPl
          nComp = (iMltPl+1)*(iMltPl+2)/2
          write(Lu,1) (Work(iAtBoMltPlAd(iMltPl)+nCenters*(iComp-1)+i*(i-1)/2+j-1),iComp=1,nComp)
        end do
        do iMltPl=0,2
          nComp = (iMltPl+1)*(iMltPl+2)/2
          write(Lu,1) (Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*(iComp-1)+i*(i-1)/2+j-1),iComp=1,nComp)
        end do
        ! Begin EB
        ! if "pointer" iAtBoPolAd is not associated with memory
        ! print zero
        if (iAtBoPolAd == 0) then
          write(Lu,1) Zero,Zero,Zero,Zero,Zero,Zero
        else
        ! End EB
          write(Lu,1) (Work(iAtBoPolAd+nCenters*k+i*(i-1)/2+j-1),k=0,5)
        ! Begin EB
        end if
        ! End EB
      end if
    end do
  end do

end if

write(Lu,6) '* Molecule properties '
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  write(Lu,1) (Work(iAtMltPlTotAd(iMltPl)+iComp-1),iComp=1,nComp)
end do
write(Lu,1) (MolPol(i),i=1,6)
write(Lu,'(A)') '* Orbital information'
write(Lu,3) NORBI,NOCOB
write(Lu,1) EneV
if (NOCOB /= 0) write(Lu,1) (OENE(I),I=1,NOCOB)
if (Method == 'UHF-SCF') then
  write(Lu,3) NOCOB_b
  if (NOcOb_b /= 0) write(Lu,1) (OENE_b(I),I=1,NOCOB_b)
end if

close(Lu)

return

1 format(3F20.10)
2 format(I5)
3 format(2I5,4X,2A)
!EB 4 format(I5,A1)
5 format(6I5,4X,A)
6 format(A)

end subroutine Wr_Files
