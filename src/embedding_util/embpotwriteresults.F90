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
! Copyright (C) 2024,2026, Lukas Schreder                              *
!***********************************************************************

subroutine embPotWriteResults(resultPath,eCASSCF,eEmb,nElec)
!***********************************************************************
!                                                                      *
! Object: routine to write the results of an embedded RASSCF run for   *
!         the FDE driver to consume.                                   *
!                                                                      *
! Called from: RASSCF                                                  *
!                                                                      *
! Author: Lukas Schreder                                               *
!                                                                      *
!***********************************************************************

  use stdalloc,    only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  implicit none

  character(len=*), intent(in) :: resultPath
  real(kind=wp), intent(in)    :: eCASSCF, eEmb, nElec

  integer(kind=iwp)           :: iunit, iatom, nAtoms
  integer(kind=iwp), external :: isFreeUnit
  real(kind=wp)               :: potNuc
  real(kind=wp), allocatable  :: coordA(:,:), chargeA(:)

  call Get_dScalar('PotNuc',potNuc)
  call Get_iScalar('Unique atoms',nAtoms)
  call mma_allocate(coordA,3,nAtoms,label='coordA')
  call mma_allocate(chargeA,nAtoms,label='chargeA')
  call Get_dArray('Unique Coordinates',coordA,3*nAtoms)
  call Get_dArray('Effective nuclear Charge',chargeA,nAtoms)

  iunit = isFreeUnit(11)
  call molcas_open(iunit,resultPath)
  write(iunit,'(ES24.14)') eCASSCF
  write(iunit,'(ES24.14)') eEmb
  write(iunit,'(ES24.14)') potNuc
  write(iunit,'(ES24.14)') nElec
  write(iunit,'(I10)') nAtoms
  do iatom=1,nAtoms
    write(iunit,'(4ES24.14)') chargeA(iatom),coordA(1,iatom),coordA(2,iatom),coordA(3,iatom)
  end do
  close(iunit)

  call mma_deallocate(coordA)
  call mma_deallocate(chargeA)

end subroutine embPotWriteResults
