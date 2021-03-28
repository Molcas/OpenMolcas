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

implicit real*8(a-h,o-z)

#include "MpData.fh"
#include "WrkSpc.fh"
#include "MolProp.fh"

parameter(mxComp=(mxMltPl+1)*(mxMltPl+2)/2)
character*16 MltPlLab, MltPlLabs(0:mxMltpl,mxComp) ! "0:" added by EB
character*16 String(0:16), PolString
real*8 MolPol(6)

iStdOut = 6

do j=1,6
  MolPol(j) = 0.0d0
end do
String(0) = 'Charge          '
String(1) = 'Dipole          '
String(2) = 'Quadrupole      '
String(3) = 'Octupole        '
String(4) = 'Hexadecapole    '
!do i=5,16
!  write(String(i),'(I2,A)') i,'-th       Cartesian'
!end do
String(5) = '5-th       Carte' !sian
String(6) = '6-th       Carte' !sian
String(7) = '7-th       Carte' !sian
String(8) = '8-th       Carte' !sian
String(9) = '9-th       Carte' !sian
String(10) = '10-th       Car' !tesian
String(11) = '11-th       Car' !tesian
String(12) = '12-th       Car' !tesian
String(13) = '13-th       Car' !tesian
String(14) = '14-th       Car' !tesian
String(15) = '15-th       Car' !tesian
String(16) = '16-th       Car' !tesian
PolString = 'Polarizability  '

if (16 < mxMltPl) then
  write(6,*) 'Increase length of MltPlLab'
  call Abend()
end if
do iMltPl=0,nMltPl
  ilab = 0
  do ix=iMltpl,0,-1
    do iy=iMltpl-ix,0,-1
      iz = iMltpl-ix-iy
      ilab = ilab+1
      do i=1,16
        MltPlLab(i:i) = ' '
      end do
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

write(iStdOut,'(1x,a16,3a16)') 'Coord           ',(MltPlLabs(1,iComp)(1:1),iComp=1,3)
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(6,'(1x,a,6a16)') String(iMltPl),(MltPlLabs(iMltPl,jComp)(1:iMltPl),jComp=iComp,min(iComp+5,nComp))
  end do
end do
if (iPol > 0) then
  write(iStdOut,'(1x,a,6a16)') PolString,(MltPlLabs(2,iComp)(1:2),iComp=1,6)
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
  write(iStdOut,'(1x,a16,3f16.8)') 'Coord           ',Cor(1,I,I),Cor(2,I,I),Cor(3,I,I)
  do iMltPl=0,nMltPl
    nComp = (iMltPl+1)*(iMltPl+2)/2
    do iComp=1,nComp,6
      write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl), &
                                       (Work(iAtBoMltPlAd(iMltPl)+nCenters*(jComp-1)+i*(i+1)/2-1),jComp=iComp,min(iComp+5,nComp))
    end do
  end do
  if (iPol > 0) then
    write(iStdOut,'(1x,a16,6f16.8)') PolString,(Work(iAtBoPolAd+nCenters*(iComp-1)+i*(i+1)/2-1),iComp=1,6)
  end if
end do
do i=1,nAtoms
  do j=1,i-1
    if (BondMat(i,j)) then
      iCount = iCount+1
      write(iStdOut,*)
      write(iStdOut,'(I5,A8,I5,A3,A10)') iCount,' Center ',i*(i-1)/2+j,'   ',CEN_LAB(i*(i-1)/2+j)
      write(iStdOut,*) '**************************************************'
      write(iStdOut,'(1x,a16,3f16.8)') 'Coord           ',Cor(1,I,j),Cor(2,I,j),Cor(3,I,j)
      do iMltPl=0,nMltPl
        nComp = (iMltPl+1)*(iMltPl+2)/2
        do iComp=1,nComp,6
          write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl), &
                                           (Work(iAtBoMltPlAd(iMltPl)+nCenters*(jComp-1)+i*(i-1)/2+j-1), &
                                            jComp=iComp,min(iComp+5,nComp))
        end do
      end do
      if (iPol > 0) then
        write(iStdOut,'(1x,a16,6f16.8)') PolString,(Work(iAtBoPolAd+nCenters*(iComp-1)+i*(i-1)/2+j-1),iComp=1,6)
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
  write(iStdOut,'(1x,a16,3f16.8)') 'Coord           ',Cor(1,I,I),Cor(2,I,I),Cor(3,I,I)
  do iMltPl=0,nMltPl
    nComp = (iMltPl+1)*(iMltPl+2)/2
    do iComp=1,nComp,6
      write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(Work(iAtMltPlAd(iMltPl)+nAtoms*(jComp-1)+i-1),jComp=iComp,min(iComp+5,nComp))
    end do
  end do
  if (iPol > 0) then
    write(iStdOut,'(1x,a16,6f16.8)') PolString,(Work(iAtPolAd+nAtoms*(iComp-1)+i-1),iComp=1,6)
    do j=1,6
      MolPol(j) = MolPol(j)+Work(iAtPolAd+nAtoms*(j-1)+i-1)
    end do
  end if
end do
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,*) ' SUMMED MULTIPOLES AND POLARIZABILITY FOR THE MOLECULE'
write(iStdOut,*) ' *****************************************************'
write(iStdOut,*)
write(iStdOut,'(1x,a16,3f16.8)') 'Coord              ',0.0d0,0.0d0,0.0d0
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(Work(iAtMltPlTotAd(iMltPl)+jComp-1),jComp=iComp,min(iComp+5,nComp))
  end do
end do
if (iPol > 0) write(iStdOut,'(1x,a16,6f16.8)') PolString,(MolPol(j),j=1,6)
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,'(1x,a16,3f16.8)') 'Coord              ',0.0d0,0.0d0,0.0d0
do iMltPl=0,nMltPl
  nComp = (iMltPl+1)*(iMltPl+2)/2
  do iComp=1,nComp,6
    write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),(Work(iAtBoMltPlTotAd(iMltPl)+jComp-1),jComp=iComp,min(iComp+5,nComp))
  end do
end do

return

!EB 10000 format(F16.6)
!EB 10003 format(3F16.6 )

end subroutine wr_mpprop
