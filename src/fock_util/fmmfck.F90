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
! Copyright (C) 2007, Mark A. Watson                                   *
!***********************************************************************

subroutine FMMFck(Dens,TwoHam,ndim)
!***********************************************************************
!                                                                      *
!     purpose: Generate FMM interface file and call FMM driver         *
!              to update Fock matrix with multipole-derived            *
!              Coulomb matrix elements                                 *
!                                                                      *
!     called from: Drv2El_dscf                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by: Mark A. Watson (maw)                                 *
!     University of Tokyo, 2007                                        *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp
#define _NOT_ACTIVE_
#ifdef _NOT_ACTIVE_
use Index_Functions, only: nTri_Elem1
use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ndim
real(kind=wp), intent(in) :: Dens(ndim), TwoHam(ndim)
#ifdef _NOT_ACTIVE_
integer(kind=iwp), parameter :: LMAX = 12
integer(kind=iwp) :: f_iostat, I, iCmp, iComp, ij, iM, iOpt, iRc, iSyLbl, iSym, J, L, lDens, M, nBas(8), nBasTot, nSym
real(kind=wp), allocatable :: CarMoms(:,:,:), CntrX(:), CntrY(:), CntrZ(:), Moms_batch(:), SphMoms(:,:,:)
character(len=8) :: Label
logical(kind=iwp) :: is_error

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)

! Compute lengths of matrices
lDens = 0
nBasTot = 0
do iSym=1,nSym
  lDens = lDens+nBas(iSym)*(nBas(iSym)+1)/2
  nBasTot = nBasTot+nBas(iSym)
end do
if (lDens /= ndim) then
  write(u6,*) 'ERROR in FMMFck',lDens,ndim
  call Abend()
end if

call mma_allocate(CntrX,ndim+4,label='CntrX')
call mma_allocate(CntrY,ndim+4,label='CntrY')
call mma_allocate(CntrZ,ndim+4,label='CntrZ')
!CntrX(:) = Zero
!CntrY(:) = Zero
!CntrZ(:) = Zero

! Read centres

iRc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'FMMCnX'
call RdOne(iRc,iOpt,Label,iComp,CntrX,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'FMMFck: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'FMMCnY'
call RdOne(iRc,iOpt,Label,iComp,CntrY,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'FMMFck: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'FMMCnZ'
call RdOne(iRc,iOpt,Label,iComp,CntrZ,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'FMMFck: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

! Read moments from one-electron files

call mma_allocate(CarMoms,ndim,nTri_Elem1(LMAX),LMAX+1,label='CarMoms')
call mma_allocate(Moms_batch,ndim+4,label='Moms_batch')
do L=0,LMAX
  do iComp=1,(L+1)*(L+2)/2
    iCmp = iComp
    iRc = -1
    iOpt = ibset(0,sNoOri)
    iSyLbl = 1
    write(Label,'(A,I2)') 'FMMInt',L
    call RdOne(iRc,iOpt,Label,iCmp,Moms_batch,iSyLbl)
    if (iRc /= 0) then
      write(u6,*) 'FMMFck: Error readin ONEINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if
    do ij=1,ndim
      CarMoms(ij,iComp,L+1) = Moms_batch(ij)
    end do
  end do
end do
call mma_deallocate(Moms_batch)

! Transform cartesian to spherical components

call mma_allocate(SphMoms,ndim,2*LMAX+1,LMAX+1,label='SphMoms')
!call fmm_call_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
call mma_deallocate(CarMoms)

! Write to FMM interface file

! Write array lengths in header file
call molcas_open_ext2(98,'MM_DATA_HEADER','SEQUENTIAL','UNFORMATTED',f_iostat,.false.,1,'REPLACE',is_error)
write(98) LMAX,nBasTot,ndim,0
close(98)

! Write multipole moments and density information
call molcas_open_ext2(98,'MM_DATA','SEQUENTIAL','UNFORMATTED',f_iostat,.false.,1,'REPLACE',is_error)

ij = 0
do J=1,nBasTot
  do I=1,J
    ij = ij+1
    do L=0,LMAX
      do M=-L,L
        iM = M+L+1
        !write (u6,'(5I3,2X,3F10.6,2E15.4)') L,M,ij,1,ij,CntrX(ij),CntrY(ij),CntrZ(ij),SphMoms(ij,iM,L+1),Dens(ij)
        write(98) L,M,I,J,ij,CntrX(ij),CntrY(ij),CntrZ(ij),SphMoms(ij,iM,L+1),Dens(ij)
      end do
    end do
  end do
end do

call mma_deallocate(CntrX)
call mma_deallocate(CntrY)
call mma_deallocate(CntrZ)
call mma_deallocate(SphMoms)

! Mark end of file with negative angular momentum
write(98) -1,0,0,0,0,Zero,Zero,Zero,Zero,Zero
close(98)

! Now call multipole code to update the Fock matrix with the
! long-range multipole-computed Coulomb matrix elements.

#include "macros.fh"
unused_var(TwoHam)
!call fmm_call_get_J_matrix(TwoHam,ndim,nBasTot,LMAX)

! Coulomb contributions of TwoHam should now be complete!

#else
#include "macros.fh"
call Untested('FMMFck')
unused_var(Dens)
unused_var(TwoHam)
#endif

return

end subroutine FMMFck
