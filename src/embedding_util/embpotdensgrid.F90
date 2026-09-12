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

subroutine embPotDensGrid(CMO,Occ,nMOs,nCMO)
!***********************************************************************
!                                                                      *
! Object: routine to write the density of RASSCF orbitals to a grid.   *
!                                                                      *
! Called from: RASSCF                                                  *
!                                                                      *
! Author: Lukas Schreder                                               *
!                                                                      *
!***********************************************************************

  use Embedding_Global, only: nEmbGridPoints, embGridCoord, embOutDensPath
  use stdalloc,         only: mma_allocate, mma_deallocate
  use Definitions,      only: wp, iwp

  implicit none

  integer(kind=iwp), intent(in) :: nMOs, nCMO
  real(kind=wp), intent(in)     :: CMO(nCMO), Occ(nMOs)

  integer(kind=iwp)              :: i, iunit, nDrv
  integer(kind=iwp), external    :: isFreeUnit
  integer(kind=iwp), allocatable :: DoIt(:)
  real(kind=wp), allocatable     :: MOValue(:), rhoGrid(:)

  ! Load grid coords: embGridCoord / nEmbGridPoints
  call EmbPotInit(.true.)

  ! initialize Seward
  i = 0
  call inisewm('mltpl',i)

  call mma_allocate(DoIt,nMOs,label='DoIt')
  DoIt(:) = 1
  call mma_allocate(MOValue,nEmbGridPoints*nMOs,label='MOval')
  call mma_allocate(rhoGrid,nEmbGridPoints,label='rhoA')

  ! evaluate density
  nDrv = 0
  call MOEval(MOValue,nMOs,nEmbGridPoints,embGridCoord,CMO,nCMO,DoIt,nDrv,1)
  call outmo(0,2,MOValue,Occ,rhoGrid,nEmbGridPoints,nMOs)

  ! write density
  iunit = isFreeUnit(11)
  call molcas_open(iunit,embOutDensPath)
  write(iunit,'(I10)') nEmbGridPoints
  do i=1,nEmbGridPoints
    write(iunit,'(ES24.14)') rhoGrid(i)
  end do
  close(iunit)

  ! cleanup
  call mma_deallocate(DoIt)
  call mma_deallocate(MOValue)
  call mma_deallocate(rhoGrid)

  call embpotfreemem()

end subroutine embPotDensGrid
