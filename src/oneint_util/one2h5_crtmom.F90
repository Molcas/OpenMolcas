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

#include "compiler_features.h"
#ifdef _HDF5_

subroutine one2h5_crtmom(fileid,nsym,nbas)
! SVC: read cartesian moments from the 1-electron integral file
! and write it to the HDF5 file specified with fileid.
! This routine does nothing if HDF5 is not supported.
! FP: also include the origins used for the operators
!
! Datasets:
!   MLTPL_X, MLTPL_Y, MLTPL_Z
!   MLTPL_XX, MLTPL_YY, MLTPL_ZZ, MLTPL_XY, MLTPL_YZ, MLTPL_XZ
!   MLTPL_ORIG

use Symmetry_Info, only: Mul
use mh5, only: mh5_init_attr, mh5_create_dset_real, mh5_put_dset,mh5_close_dset

implicit none
integer :: fileid
integer :: nsym, nbas(*)
#include "Molcas.fh"
#include "stdalloc.fh"
integer :: isym, jsym, msym
integer :: nb, nbast, nB1, nB2
real*8, allocatable :: MLTPL(:,:), Scratch(:)
real*8, dimension(3,3) :: mp_orig
integer :: iRc, iOpt, iComp, iSyMsk
character(len=8) :: Label
character(len=1) :: mltpl1_comp(3) = ['X','Y','Z']
character(len=2) :: mltpl2_comp(6) = ['XX','XY','XZ','YY','YZ','ZZ']
integer :: i, j, iOff, jOff, iScrOff, iBas, jBas
integer :: dsetid

nbast = 0
do isym=1,nsym
  nb = nbas(isym)
  nbast = nbast+nb
end do

mp_orig(:,:) = 0.

call mma_allocate(MLTPL,NBAST,NBAST)
call mma_allocate(Scratch,NBAST**2+3)

do icomp=1,3
  MLTPL = 0.0d0
  iRc = -1
  iOpt = 4
  iSyMsk = 0
  Label = 'Mltpl  1'
  call RdOne(iRc,iOpt,Label,iComp,Scratch,iSyMsk)
  ! iSyMsk tells us which symmetry combination is valid
  iScrOff = 0
  iOff = 0
  do iSym=1,nSym
    jOff = 0
    nB1 = nBas(iSym)
    do jSym=1,iSym
      mSym = Mul(iSym,jSym)
      nB2 = nBas(jSym)
      if (iand(2**(mSym-1),iSyMsk) /= 0) then
        if (iSym == jSym) then
          do j=1,nB2
            jBas = jOff+j
            do i=1,j
              iBas = iOff+i
              MLTPL(iBas,jBas) = Scratch(1+iScrOff)
              iScrOff = iScrOff+1
            end do
          end do
        else
          do j=1,nB2
            jBas = jOff+j
            do i=1,nB1
              iBas = iOff+i
              MLTPL(jBas,iBas) = Scratch(1+iScrOff)
              iScrOff = iScrOff+1
            end do
          end do
        end if
      end if
      jOff = jOff+nB2
    end do
    iOff = iOff+nB1
  end do
  do j=1,nBasT
    do i=1,j-1
      MLTPL(j,i) = MLTPL(i,j)
    end do
  end do
  dsetid = mh5_create_dset_real(fileid,'AO_MLTPL_'//mltpl1_comp(icomp),2,[NBAST,NBAST])
  call mh5_init_attr(dsetid,'DESCRIPTION', &
                     '1st-order multipole matrix of the atomic orbitals, arranged as matrix of size [NBAST,NBAST]')
  call mh5_put_dset(dsetid,MLTPL)
  call mh5_close_dset(dsetid)
end do

mp_orig(1:3,2) = Scratch(iScrOff+1:iScrOff+3)

do icomp=1,6
  MLTPL = 0.0d0
  iRc = -1
  iOpt = 4
  iSyMsk = 0
  Label = 'Mltpl  2'
  call RdOne(iRc,iOpt,Label,iComp,Scratch,iSyMsk)
  ! iSyMsk tells us which symmetry combination is valid
  iScrOff = 0
  iOff = 0
  do iSym=1,nSym
    jOff = 0
    nB1 = nBas(iSym)
    do jSym=1,iSym
      mSym = Mul(iSym,jSym)
      nB2 = nBas(jSym)
      if (iand(2**(mSym-1),iSyMsk) /= 0) then
        if (iSym == jSym) then
          do j=1,nB2
            jBas = jOff+j
            do i=1,j
              iBas = iOff+i
              MLTPL(iBas,jBas) = Scratch(1+iScrOff)
              iScrOff = iScrOff+1
            end do
          end do
        else
          do j=1,nB2
            jBas = jOff+j
            do i=1,nB1
              iBas = iOff+i
              MLTPL(jBas,iBas) = Scratch(1+iScrOff)
              iScrOff = iScrOff+1
            end do
          end do
        end if
      end if
      jOff = jOff+nB2
    end do
    iOff = iOff+nB1
  end do
  do j=1,nBasT
    do i=1,j-1
      MLTPL(j,i) = MLTPL(i,j)
    end do
  end do
  dsetid = mh5_create_dset_real(fileid,'AO_MLTPL_'//mltpl2_comp(icomp),2,[NBAST,NBAST])
  call mh5_init_attr(dsetid,'DESCRIPTION', &
                     '2nd-order multipole matrix of the atomic orbitals, arranged as matrix of size [NBAST,NBAST]')
  call mh5_put_dset(dsetid,MLTPL)
  call mh5_close_dset(dsetid)
end do

mp_orig(1:3,3) = Scratch(iScrOff+1:iScrOff+3)

call mma_deallocate(MLTPL)
call mma_deallocate(Scratch)

dsetid = mh5_create_dset_real(fileid,'MLTPL_ORIG',2,[3,3])
call mh5_init_attr(dsetid,'DESCRIPTION','Origin used for the multipole moment operators: arranged as overlap, dipole, quadrupole')
call mh5_put_dset(dsetid,mp_orig,[3,3],[0,0])
call mh5_close_dset(dsetid)

end subroutine one2h5_crtmom

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(one2h5_crtmom)

#endif
