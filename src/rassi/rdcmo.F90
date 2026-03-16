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

subroutine RDCMO_RASSI(JOB,CMO)

#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_dset, mh5_close_file
#endif
use rassi_aux, only: ipglob
use stdalloc, only: mma_allocate, mma_deallocate
use Cntrl, only: NJOB, PRORB, JBNAME
use cntrl, only: IDCMO, iTOC15, LuIph
use Symmetry_Info, only: nSym => nIrrep
use rassi_data, only: NCMO, NBASF, NOSH
use Constants, only: Zero
use Definitions, only: u6

implicit none
#ifdef _HDF5_
integer :: refwfn_id
#endif
integer JOB
real*8 CMO(NCMO)
integer I, IAD, IDISK, ISY, L1, L2, LEN, NB, NBUF
real*8, allocatable :: BUF(:)

CMO(:) = Zero
if ((JOB < 1) .or. (JOB > NJOB)) then
  write(u6,*) ' RDCMO_RASSI: Invalid JOB parameter.'
  write(u6,*) ' JOB, NJOB:',JOB,NJOB
  call ABEND()
end if
if (IPGLOB >= 4) write(u6,*) ' RDCMO_RASSI called for file '//trim(JBNAME(JOB))
! READ ORBITAL COEFFICIENTS FROM INTERFACE. ORIGINALLY ALL
! CMO COEFFS, INCLUDING VIRTUALS, WERE WRITTEN CONTIGUOUSLY.
NBUF = 0
do I=1,NSYM
  NBUF = NBUF+NBASF(I)**2
end do
call mma_allocate(BUF,NBUF,Label='BUF')

#ifdef _HDF5_
if (mh5_is_hdf5(jbname(job))) then
  !*********************************************************************
  !
  ! For HDF5 formatted job files
  !
  !*********************************************************************
  refwfn_id = mh5_open_file_r(jbname(job))
  call mh5_fetch_dset(refwfn_id,'MO_VECTORS',BUF)
  call mh5_close_file(refwfn_id)
else
#endif
  !*********************************************************************
  !
  ! For JOBIPH/JOBMIX formatted job files
  !
  !*********************************************************************
  call DANAME(LUIPH,JBNAME(JOB))
  IAD = 0
  call IDAFILE(LUIPH,2,ITOC15,30,IAD)
  IDISK = IDCMO(JOB)
  call DDAFILE(LUIPH,2,BUF,NBUF,IDISK)
  call DACLOS(LUIPH)
#ifdef _HDF5_
end if
!***********************************************************************
#endif

if (IPGLOB >= 5) then
  write(u6,*) ' Reading CMO'
  write(u6,*) ' NBUF=',NBUF
  write(u6,*) ' Array read in:'
  write(u6,'(1x,5f16.8)') (buf(i),i=1,nbuf)
end if
L1 = 1
L2 = 1
do ISY=1,NSYM
  NB = NBASF(ISY)
  LEN = NOSH(ISY)*NB
  if (LEN > 0) call DCOPY_(LEN,BUF(L1),1,CMO(L2),1)
  L2 = L2+LEN
  L1 = L1+NB**2
end do
if (IPGLOB >= 5) then
  write(u6,*) ' Gathered CMO from array.'
  write(u6,*) ' NCMO=',NCMO
  write(u6,'(1x,5f16.8)') (CMO(i),i=1,ncmo)
end if
call mma_deallocate(BUF)

if ((IPGLOB > 0) .and. PRORB) then
  write(u6,*)
  call WRMAT('MO ORBITAL COEFFICIENTS:',1,NBASF,NOSH,NCMO,CMO)
end if

end subroutine RDCMO_RASSI
