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
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  RdVec_HDF5
!
!> @brief
!>   Read orbital data from an HDF5 file
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Similar to ::RdVec, but for HDF5 files. The \p Label argument can
!>  contain a `B` if one wants the data for beta orbitals.
!>
!> @param[in]  fileid   Identifier of the open HDF5 file
!> @param[in]  Label    Properties to read from the file
!> @param[in]  nSym     Number of irreps
!> @param[in]  nBas     Number of basis functions per irrep
!> @param[out] CMO      Orbital coefficients
!> @param[out] Occ      Orbital occupations
!> @param[out] Ene      Orbital energies
!> @param[out] Ind      Orbital type indices
!***********************************************************************

subroutine RdVec_HDF5(fileid,Label,nSym,nBas,CMO,Occ,Ene,Ind)

#ifdef _HDF5_
use mh5, only: mh5_exists_dset, mh5_fetch_dset
#endif

implicit none
integer, intent(In) :: fileid, nSym, nBas(nSym)
character(Len=*), intent(In) :: Label
real*8, dimension(*) :: CMO, Occ, Ene
integer, dimension(*) :: Ind
#ifdef _HDF5_
#include "stdalloc.fh"
integer :: Beta, nB
character(Len=128) :: DataSet, su, sl
character(Len=1), allocatable :: typestring(:)

Beta = 0
su = ''
sl = ''
if (index(Label,'A') > 0) then
  Beta = -1
  su = 'ALPHA_'
  sl = 'alpha '
end if
if (index(Label,'B') > 0) then
  if (Beta /= 0) then
    write(6,*)
    call AbEnd()
  end if
  Beta = 1
  su = 'BETA_'
  sl = 'beta '
end if

if (index(Label,'E') > 0) then
  DataSet = 'MO_'//trim(su)//'ENERGIES'
  if (mh5_exists_dset(fileid,DataSet)) then
    call mh5_fetch_dset(fileid,DataSet,Ene)
  else
    write(6,*) 'The HDF5 file does not contain '//trim(sl)//'MO energies.'
    call AbEnd()
  end if
end if

if (index(Label,'O') > 0) then
  DataSet = 'MO_'//trim(su)//'OCCUPATIONS'
  if (mh5_exists_dset(fileid,DataSet)) then
    call mh5_fetch_dset(fileid,DataSet,Occ)
  else
    write(6,*) 'The HDF5 file does not contain '//trim(sl)//'MO occupations.'
    call AbEnd()
  end if
end if

if (index(Label,'C') > 0) then
  DataSet = 'MO_'//trim(su)//'VECTORS'
  if (mh5_exists_dset(fileid,DataSet)) then
    call mh5_fetch_dset(fileid,DataSet,CMO)
  else
    write(6,*) 'The HDF5 file does not contain '//trim(sl)//'MO coefficients.'
    call AbEnd()
  end if
end if

if (index(Label,'I') > 0) then
  nB = sum(nBas)
  call mma_allocate(typestring,nB)
  DataSet = 'MO_'//trim(su)//'TYPEINDICES'
  if (mh5_exists_dset(fileid,DataSet)) then
    call mh5_fetch_dset(fileid,DataSet,typestring)
    call tpstr2tpidx(typestring,Ind,nB)
  end if
  call mma_deallocate(typestring)
end if

#else
call WarningMessage(2,'Calling RdVec_HDF5, but HDF5 is disabled')
call AbEnd()

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(fileid)
  call Unused_character(Label)
  call Unused_integer(nSym)
  call Unused_integer_array(nBas)
  call Unused_real_array(CMO)
  call Unused_real_array(Occ)
  call Unused_real_array(Ene)
  call Unused_integer_array(Ind)
end if
#endif

end subroutine RdVec_HDF5
