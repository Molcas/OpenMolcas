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
! Copyright (C) 2002, Laura Gagliardi                                  *
!               2002, Roland Lindh                                     *
!***********************************************************************

subroutine LoProp(ireturn)
!***********************************************************************
!      Author:Laura Gagliardi, Dipartimento di Chimica Fisica,         *
!             University of Palermo, ITALY. December 2002              *
!             Roland Lindh, Department of Chemical Physics,            *
!             University of Lund, SWEDEN.                              *
!***********************************************************************

use loprop_arrays, only: LP_context_type
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: ireturn
#include "itmax.fh"
#include "Molcas.fh"
integer(kind=iwp), parameter :: nElem=(iTabMx+1)*(iTabMx+2)*(iTabMx+3)/6
integer(kind=iwp) :: i, iPert, iPL, iPlot, iPrint, lMax, mElem, nAtoms, nBas(8), nBas1, nBas2, nBasMax, nij, nOrb(8), nPert, &
                     nSize, nStateF, nStateI, nSym, nTemp
real(kind=wp) :: Bond_Threshold, CoC(3), Origin(3,0:iTabMx), SubScale
logical(kind=iwp) :: lSave, NoField, PrintDen, Restart, Standard, SubtractDen, TDensity, UserDen, Utility
character(len=LenIn4) :: LblCnt(MxAtom)
character(len=12) :: Opt_Method
type(LP_context_type) :: LP_context
type(Alloc1DArray_Type) :: D(0:6), imu(0:nElem-1)
real(kind=wp), allocatable :: EC(:,:), MP(:,:,:), MPq(:), tmp(:), sq_mu(:,:), sq_temp(:), Ttot(:,:), Ttot_Inv(:,:)
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
Utility = .false.
lSave = ireturn == 0
ireturn = 99

iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

!                                                                      *
!***********************************************************************
!                                                                      *
! Prelude

call Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,LP_context,nSize,nBas1,nBas2,nBasMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the input
!
! NoField is defaulted to true if symmetry is implied.

NoField = .true.
Standard = .true.
UserDen = .false.
PrintDen = .false.
SubtractDen = .false.
Restart = .false.
TDensity = .false.
nStateI = 1
nStateF = 1
Bond_Threshold = -One
iPlot = 0
iPrint = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Do the LoProp localization.

call mma_allocate(Ttot,nBas1,nBas1,label='Ttot')
call mma_allocate(Ttot_Inv,nBas1,nBas1,label='TtotInv')

call Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,LP_context%center,LP_context%otype,nBas1,nBas2,nSym,nBasMax,LP_context%PInv,Restart)

call mma_deallocate(LP_context%otype)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read in the multipole moment integrals once and for all.

lMax = 0   ! do only charges
mElem = (lMax*(lMax**2+6*lMax+11)+6)/6

nTemp = nBas1**2
call mma_allocate(tmp,nTemp,label='tmp')

call mma_allocate(sq_mu,[1,nBas1**2],[0,mElem-1],label='sq_mu')
call mma_allocate(MPq,mElem,label='MPq')
call Read_Multipole_Int(lMax,sq_mu,nBas,imu,Ttot,tmp,Origin,MPq,mElem,nBas1,nBas2,nBasMax,nTemp,nSym,LP_context%PInv,Restart, &
                        Utility)
call mma_deallocate(Ttot)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the 1-particle density matrix

iPert = 0
call Get_Density_Matrix(D(0),nBas1,nBas2,nBasMax,nBas,nSym,LP_context%P,UserDen,PrintDen,SubtractDen,SubScale,LP_context%Q_Nuc, &
                        nAtoms,iPert,Restart,Utility,TDensity,nStateI,nStateF)

do i=mElem-1,0,-1
  call mma_deallocate(imu(i)%A)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!     S T A T I C   P R O P E R T I E S                                *
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute localized multipole moments

nPert = 2*3+1
if (NoField) nPert = 1
nij = (nAtoms*(nAtoms+1)/2)
call mma_allocate(MP,nij,mElem,nPert,label='MP')
call mma_allocate(sq_temp,nTemp,label='sq_temp')
call mma_allocate(EC,3,nij,label='EC')

if (iPL >= 2) then
  write(u6,*)
  call CollapseOutput(1,'   Static properties:')
  write(u6,'(3X,A)') '   ------------------'
  write(u6,*)
end if
call Local_Properties(LP_context%C,nAtoms,sq_mu,mElem,sq_temp,Origin,LP_context%center,Ttot_Inv,tmp,nij,nPert,D,MP,lMax,MPq,CoC, &
                      EC,LP_context%ANr,Standard,nBas1,nTemp,LP_context%Q_Nuc,Bond_Threshold,Opt_Method,iPlot,iPrint,nSym)

call mma_deallocate(sq_mu)
call mma_deallocate(Ttot_Inv)
call mma_deallocate(sq_temp)
call mma_deallocate(LP_context%center)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the properties

call Get_cArray('LP_L',LblCnt,LenIn4*nAtoms)
call LoProp_Print(MP,nij,nElem,nAtoms,LP_context%Q_Nuc,LblCnt,lSave)
if (iPL >= 2) then
  call CollapseOutput(0,'   Static properties:')
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(MPq)
call mma_deallocate(EC)
call mma_deallocate(MP)
call mma_deallocate(tmp)
call mma_deallocate(LP_context%Q_Nuc)
call mma_deallocate(LP_context%ANr)
call mma_deallocate(LP_context%C)
call mma_deallocate(LP_context%P)
call mma_deallocate(LP_context%PInv)
!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = 0
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine LoProp
