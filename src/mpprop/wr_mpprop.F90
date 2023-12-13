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

subroutine wr_mpprop(nAtoms,nCenters,nMltPl,iPol)

use MPprop_globals, only: AtBoMltPl, AtBoMltPlTot, AtMltPl, AtMltPlTot, AtPol, AtBoPol, BondMat, Cen_Lab, Cor, Title
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nCenters, nMltPl, iPol
integer(kind=iwp) :: i, iComp, iCount, ilab, iMltPl, iStdOut, ix, ixx, iy, iyy, iz, izz, j, jComp, mxComp, nComp
real(kind=wp) :: MolPol(6)
character(len=16) :: MltPlLab, PolString
character(len=16), allocatable :: MltPlLabs(:,:), String(:)

iStdOut = u6

mxComp = (nMltPl+1)*(nMltpl+2)/2
call mma_allocate(MltPlLabs,[0,max(nMltPl,2)],[1,mxComp],label='MltPlLabs')
call mma_allocate(String,[0,max(nMltPl,4)],label='String')

MolPol(:) = Zero
String(0) = 'Charge'
String(1) = 'Dipole'
String(2) = 'Quadrupole'
String(3) = 'Octupole'
String(4) = 'Hexadecapole'
do i=5,nMltPl
  write(String(i),'(i4,"th")') i
  if (String(i)(3:3) /= '1') then
    select case (String(i)(4:4))
      case ('1')
        String(i)(5:6) = 'st'
      case ('2')
        String(i)(5:6) = 'nd'
      case ('3')
        String(i)(5:6) = 'rd'
      case default
    end select
  end if
  String(i) = adjustl(String(i))
  String(i)(8:16) = 'Cartesian'
end do
PolString = 'Polarizability'

do iMltPl=0,max(nMltPl,2)
  ilab = 0
  do ix=iMltpl,0,-1
    do iy=iMltpl-ix,0,-1
      iz = iMltpl-ix-iy
      ilab = ilab+1
      MltPlLab = ' '
      do izz=iMltPl,iMltPl-iz+1,-1
        MltPlLab(izz:izz) = 'Z'
      end do
      do iyy=iMltPl-iz,iMltPl-iz-iy+1,-1
        MltPlLab(iyy:iyy) = 'Y'
      end do
      !max( ,1) added by EB
      do ixx=iMltPl-iz-iy,max(iMltPl-nMltPl+1,1),-1
        MltPlLab(ixx:ixx) = 'X'
      end do
      MltPlLabs(iMltPl,ilab) = MltPlLab
    end do
  end do
end do

write(iStdOut,*)
write(iStdOut,*) ' The name of the molecule will be',Title
write(iStdOut,*)
write(iStdOut,*) ' THE MULTIPOLES AND POLARIZABILITY FOR THE EXPANSION'
write(iStdOut,*) ' ****************************************************************'
write(iStdOut,*)
write(iStdOut,*) ' TOTAL NUMBER OF CENTERS   : ',nCenters
write(iStdOut,*) ' OF WHICH THERE ARE ATOMS  : ',nAtoms
write(iStdOut,*) ' AND BONDS IN THE MOLECULE : ',nCenters-nAtoms
write(iStdOut,*)

write(iStdOut,'(1x,a16,3a16)') 'Coord',(MltPlLabs(1,iComp),iComp=1,3)
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(u6,'(1x,a,6a16)') String(iMltPl),(MltPlLabs(iMltPl,jComp),jComp=iComp,min(iComp+5,nComp))
  end do
end do
if (iPol > 0) then
  write(iStdOut,'(1x,a,6a16)') PolString,(MltPlLabs(2,iComp),iComp=1,6)
end if
iCount = 0
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*) 'Multipole expansion for Atoms and Bonds'
write(iStdOut,*) '***************************************'
do i=1,nAtoms
  iCount = iCount+1
  write(iStdOut,*)
  write(iStdOut,'(I5,A8,I5,A3,A10)') iCount,' Center ',i*(i+1)/2,'   ',CEN_LAB(i*(i+1)/2)
  write(iStdOut,*) '**************************************************'
  write(iStdOut,'(1x,a16,3f16.8)') 'Coord',Cor(1,I,I),Cor(2,I,I),Cor(3,I,I)
  do iMltPl=0,nMltPl
    nComp = (iMltPl+1)*(iMltPl+2)/2
    do iComp=1,nComp,6
      write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(AtBoMltPl(iMltPl)%A(jComp,i*(i+1)/2),jComp=iComp,min(iComp+5,nComp))
    end do
  end do
  if (iPol > 0) then
    write(iStdOut,'(1x,a16,6f16.8)') PolString,AtBoPol(:,i*(i+1)/2)
  end if
end do
do i=1,nAtoms
  do j=1,i-1
    if (BondMat(i,j)) then
      iCount = iCount+1
      write(iStdOut,*)
      write(iStdOut,'(I5,A8,I5,A3,A10)') iCount,' Center ',i*(i-1)/2+j,'   ',CEN_LAB(i*(i-1)/2+j)
      write(iStdOut,*) '**************************************************'
      write(iStdOut,'(1x,a16,3f16.8)') 'Coord',Cor(1,i,j),Cor(2,i,j),Cor(3,i,j)
      do iMltPl=0,nMltPl
        nComp = (iMltPl+1)*(iMltPl+2)/2
        do iComp=1,nComp,6
          write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(AtBoMltPl(iMltPl)%A(jComp,i*(i-1)/2+j),jComp=iComp,min(iComp+5,nComp))
        end do
      end do
      if (iPol > 0) then
        write(iStdOut,'(1x,a16,6f16.8)') PolString,AtBoPol(:,i*(i-1)/2+j)
      end if
    end if
  end do
end do
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*) 'Multipole expansion for Atoms'
write(iStdOut,*) '*****************************'

do i=1,nAtoms
  write(iStdOut,*)
  write(iStdOut,'(I5,A8,I5,A3,A10)') i,' Atom   ',i,'   ',CEN_LAB(i*(i+1)/2)
  write(iStdOut,*) '**************************************************'
  write(iStdOut,'(1x,a16,3f16.8)') 'Coord',Cor(1,i,i),Cor(2,i,i),Cor(3,i,i)
  do iMltPl=0,nMltPl
    nComp = (iMltPl+1)*(iMltPl+2)/2
    do iComp=1,nComp,6
      write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(AtMltPl(iMltPl)%A(jComp,i),jComp=iComp,min(iComp+5,nComp))
    end do
  end do
  if (iPol > 0) then
    write(iStdOut,'(1x,a16,6f16.8)') PolString,AtPol(:,i)
    MolPol(:) = MolPol(:)+AtPol(:,i)
  end if
end do
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*) ' SUMMED MULTIPOLES AND POLARIZABILITY FOR THE MOLECULE'
write(iStdOut,*) ' *****************************************************'
write(iStdOut,*)
write(iStdOut,'(1x,a16,3f16.8)') 'Coord',Zero,Zero,Zero
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(AtMltPlTot(iMltPl)%A(jComp),jComp=iComp,min(iComp+5,nComp))
  end do
end do
if (iPol > 0) write(iStdOut,'(1x,a16,6f16.8)') PolString,(MolPol(j),j=1,6)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,'(1x,a16,3f16.8)') 'Coord',Zero,Zero,Zero
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(AtBoMltPlTot(iMltPl)%A(jComp),jComp=iComp,min(iComp+5,nComp))
  end do
end do

call mma_deallocate(MltPlLabs)
call mma_deallocate(String)

return

!EB 10000 format(F16.6)
!EB 10003 format(3F16.6 )

end subroutine wr_mpprop
