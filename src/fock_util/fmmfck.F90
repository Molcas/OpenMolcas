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
#ifdef _NOT_ACTIVE_
use Constants, only: Zero
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ndim
real(kind=wp) :: Dens(ndim), TwoHam(ndim)
#ifdef _NOT_ACTIVE_
! Define local variables
!#include "mxdm.fh"
integer(kind=iwp), parameter :: LMAX = 12
integer(kind=iwp) :: I, iComp, ij, iM, iOpt, iRc, iSyLbl, iSym, J, L, lDens, M, nBas(8), nBasTot, nSym
real(kind=wp) :: CarMoms(ndim,(LMAX+1)*(LMAX+2)/2,LMAX+1), CntrX(ndim+4), CntrY(ndim+4), CntrZ(ndim+4), Moms_batch(ndim+4), &
                 SphMoms(ndim,2*LMAX+1,LMAX+1)
character(len=8) :: Label

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

!call dcopy_(ndim,Zero,0,Moments,1)
!call dcopy_(ndim,Zero,0,CntrX,1)
!call dcopy_(ndim,Zero,0,CntrY,1)
!call dcopy_(ndim,Zero,0,CntrZ,1)

! Read centres

iRc = -1
iOpt = 2
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

do L=0,LMAX
  do iComp=1,(L+1)*(L+2)/2
    iRc = -1
    iOpt = 2
    iSyLbl = 1
    write(Label,'(A,I2)') 'FMMInt',L
    call RdOne(iRc,iOpt,Label,iComp,Moms_batch,iSyLbl)
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

! Transform cartesian to spherical components

!call fmm_call_car_to_sph(CarMoms,SphMoms,ndim,LMAX)

! Write to FMM interface file

! Write array lengths in header file
open(98,FILE='MM_DATA_HEADER',FORM='UNFORMATTED',STATUS='REPLACE')
write(98) LMAX,nBasTot,ndim,0
close(98,STATUS='KEEP')

! Write multipole moments and density information
open(98,FILE='MM_DATA',FORM='UNFORMATTED',STATUS='REPLACE')

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

! Mark end of file with negative angular momentum
write(98)-1,0,0,0,0,Zero,Zero,Zero,Zero,Zero
close(98,STATUS='KEEP')

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
