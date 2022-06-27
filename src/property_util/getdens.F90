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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine GetDens(FName,Density,iPrint)
!***********************************************************************
! Object: to get the 1 particle density from file INPORB               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chem. Phys.                       *
!             University of Lund, SWEDEN                               *
!             January 2000                                             *
!***********************************************************************

use PrpPnt, only: Den, nDen, nOcc, nVec, Occ, Vec
use Basis_Info, only: nBas
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_close_file
#endif
use stdalloc, only: mma_allocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: FName
logical(kind=iwp), intent(in) :: Density
integer(kind=iwp), intent(in) :: iPrint
integer(kind=iwp) :: iadDen, iadOcc, iadVec, iBas, icd, ico, ictd, icto, ictv, icv, iDummy(1), iErr, iIrrep, j1, j2, LuVec
#ifdef _HDF5_
integer(kind=iwp) :: id_file
#endif
real(kind=wp) :: Dummy(1)
character(len=80) :: Line
integer(kind=iwp), external :: n2Tri

nDen = n2Tri(1)
nVec = S%n2Tot
nOcc = S%nDim
if (Density) call mma_allocate(Den,nDen,label='Den')
iadDen = 1
call mma_allocate(Vec,nVec,label='Vec')
call mma_allocate(Occ,nOcc,label='Occ')
iadDen = 1
iadVec = 1
iadOcc = 1
#ifdef _HDF5_
if (mh5_is_hdf5(FName)) then
  id_file = mh5_open_file_r(FName)
  call RdVec_HDF5(id_file,'CO',nIrrep,nBas,Vec,Occ,Dummy,iDummy)
  call mh5_close_file(id_file)
  write(u6,*)
  write(u6,'(A,1X,A)') ' Vectors read from HDF5 file:',trim(FName)
  write(u6,*)
else
#endif
  LuVec = 19
  call RdVec(FName,LuVec,'CO',nIrrep,nBas,nBas,Vec,Occ,Dummy,iDummy,Line,0,iErr)
  write(u6,*)
  write(u6,'(A)') ' Header from vector file:'
  write(u6,*)
  write(u6,'(A)') trim(Line)
  write(u6,*)
#ifdef _HDF5_
end if
#endif

if (Density) then

  ! Build the density matrix.

  Den(:) = Zero

  ictv = iadVec
  icto = iadOcc
  ictd = iadDen
  do iIrrep=0,nIrrep-1
    do iBas=1,nBas(iIrrep)
      icv = ictv
      ico = icto
      icd = ictd
      do j1=0,nBas(iIrrep)-1
        do j2=0,j1-1
          Den(icd) = Den(icd)+Occ(ico)*Vec(icv+j1)*Two*Vec(icv+j2)
          icd = icd+1
        end do
        Den(icd) = Den(icd)+Occ(ico)*Vec(icv+j1)*Vec(icv+j1)
        icd = icd+1
      end do
      ictv = ictv+nBas(iIrrep)
      icto = icto+1
    end do
    ictd = ictd+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end do
  iadVec = iadDen
  iadOcc = iadDen
  nOcc = nDen
  nVec = nDen
  if (iPrint >= 10) call PrMtrx(' Density matrix',[1],1,[iadDen],Den)

end if

return

end subroutine GetDens
