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

subroutine one2h5_crtmom(fileid,nSym,nBas)
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
use OneDat, only: sNoNuc
use mh5, only: mh5_close_dset, mh5_create_dset_real, mh5_init_attr, mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: fileid, nSym, nBas(nSym)
integer(kind=iwp) :: dsetid, i, iBas, iCmp, iComp, iOff, iOpt, iRc, iScrOff, iSym, iSyMsk, j, jBas, jOff, jsym, msym, nB1, nB2, &
                     nbast
real(kind=wp) :: mp_orig(3,3)
character(len=8) :: Label
real(kind=wp), allocatable :: MLTPL(:,:), Scratch(:)
character(len=*), parameter :: mltpl1_comp(3) = ['X','Y','Z'], mltpl2_comp(6) = ['XX','XY','XZ','YY','YZ','ZZ']

nbast = 0
do iSym=1,nSym
  nbast = nbast+nBas(iSym)
end do

mp_orig(:,:) = 0.

call mma_allocate(MLTPL,NBAST,NBAST)
call mma_allocate(Scratch,NBAST**2+3)

do icomp=1,3
  iCmp = iComp
  MLTPL = Zero
  iRc = -1
  iOpt = ibset(0,sNoNuc)
  iSyMsk = 0
  Label = 'Mltpl  1'
  call RdOne(iRc,iOpt,Label,iCmp,Scratch,iSyMsk)
  ! iSyMsk tells us which symmetry combination is valid
  iScrOff = 0
  iOff = 0
  do iSym=1,nSym
    jOff = 0
    nB1 = nBas(iSym)
    do jSym=1,iSym
      mSym = Mul(iSym,jSym)
      nB2 = nBas(jSym)
      if (btest(iSyMsk,mSym-1)) then
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
  iCmp = iComp
  MLTPL = Zero
  iRc = -1
  iOpt = ibset(0,sNoNuc)
  iSyMsk = 0
  Label = 'Mltpl  2'
  call RdOne(iRc,iOpt,Label,iCmp,Scratch,iSyMsk)
  ! iSyMsk tells us which symmetry combination is valid
  iScrOff = 0
  iOff = 0
  do iSym=1,nSym
    jOff = 0
    nB1 = nBas(iSym)
    do jSym=1,iSym
      mSym = Mul(iSym,jSym)
      nB2 = nBas(jSym)
      if (btest(iSyMsk,mSym-1)) then
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

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(one2h5_crtmom)

#endif
