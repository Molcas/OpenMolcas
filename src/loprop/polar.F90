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

subroutine Polar(ireturn)
!***********************************************************************
!      Author:Laura Gagliardi, Dipartimento di Chimica Fisica,         *
!             University of Palermo, ITALY. December 2002              *
!             Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN.                                         *
!***********************************************************************

use loprop_arrays, only: LP_context_type
use Data_Structures, only: Allocate_DT, Alloc1DArray_Type, Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "Molcas.fh"
integer(kind=iwp) :: i, i_f, iPert, iPlot, iPrint, lMax, LoProp_Mode, LuYou, mElem, MpProp_Level, nAtoms, nBas(8), nBas1, nBas2, &
                     nBasMax, nDim, nij, nOcOb, nOrb(8), nPert, nSize, nStateF, nStateI, nSym, nTemp, nThrs
real(kind=wp) :: Alpha, Bond_Threshold, CoC(3), Delta, Dlt, dLimmo(2), Energy_Without_FFPT, Ep, SubScale, Thrs1, Thrs2, ThrsMul
logical(kind=iwp) :: NoField, Standard, Utility, UserDen, PrintDen, SubtractDen, Restart, TDensity, Diffuse(3), Exists, LIonize
character(len=12) :: Opt_Method
type(LP_context_type) :: LP_context
type(Alloc1DArray_Type) :: D(0:6)
character(len=LenIn), allocatable :: LblCnt(:)
character(len=LenIn4), allocatable :: LblCnt4(:)
real(kind=wp), allocatable :: Cpl(:,:), CplT(:,:), EC(:,:), Ene_Occ(:), h0(:), MP(:,:,:), MPp(:,:), MPq(:), nxMP(:), Origin(:,:), &
                              Pol(:,:), sq_mu(:,:), sq_temp(:), tmp(:), TP(:), Ttot(:,:), Ttot_Inv(:,:), xMP(:), xxMP(:)
type(Alloc1DArray_Type), allocatable :: imu(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
Utility = .false.
ireturn = 99
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

NoField = nSym /= 1
call Readin_polar(NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen,SubtractDen, &
                  SubScale,Restart,TDensity,nStateI,nStateF,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,LIonize)
call mma_allocate(Ene_Occ,nBas1,label='Ene_Occ')
call InfoToMp(nSym,nBas,nBas1,Energy_Without_FFPT,Ene_Occ,nOcOb,UserDen,Restart)
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

call Get_iScalar('Highest Mltpl',lMax)
write(u6,'(A,I2)') ' Multipole moments will be processed up to order ',lMax
write(u6,*)
mElem = (lMax+1)*(lMax+2)*(lMax+3)/6

nTemp = nBas1**2
call mma_allocate(tmp,nTemp,label='tmp')

call mma_allocate(Origin,[1,3],[0,lMax],label='Origin')
call mma_allocate(sq_mu,[1,nBas1**2],[0,mElem-1],label='sq_mu')
call mma_allocate(MPq,mElem,label='MPq')
call Allocate_DT(imu,[0,mElem-1],label='imu')
call Read_Multipole_Int(lMax,sq_mu,nBas,imu,Ttot,tmp,Origin,MPq,mElem,nBas1,nBas2,nBasMax,nTemp,nSym,LP_context%PInv,Restart, &
                        Utility)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the 1-particle density matrix

Dlt = -Delta
iPert = 0
call Get_Density_Matrix(D(0),nBas1,nBas2,nBasMax,nBas,nSym,LP_context%P,UserDen,PrintDen,SubtractDen,SubScale,LP_context%Q_Nuc, &
                        nAtoms,iPert,Restart,Utility,TDensity,nStateI,nStateF)

if (.not. NoField) then
  ! This assumes mElem >= 4, so lMax >= 1
  if (lMax < 1) call Abend()
  ! Read the one-electron hamiltonian.
  call mma_allocate(h0,nSize,label='h0')
  call Read_h0(nSize,h0,Restart)
  do iPert=1,6
    i_f = (iPert+1)/2
    Dlt = -Dlt
    if ((.not. Restart) .and. (.not. UserDen)) then
      call Comp_F(h0,imu(i_f)%A,nBas(1),Dlt,Ep,imu(0)%A,CoC(i_f),Origin(i_f,1))
    end if
    call Get_Density_Matrix(D(iPert),nBas1,nBas2,nBasMax,nBas,nSym,LP_context%P,UserDen,PrintDen,SubtractDen,SubScale, &
                            LP_context%Q_Nuc,nAtoms,iPert,Restart,Utility,TDensity,nStateI,nStateF)
  end do
  call mma_deallocate(h0)
end if
call Deallocate_DT(imu)
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

call Local_Properties(LP_context%C,nAtoms,sq_mu,mElem,sq_temp,Origin,LP_context%center,Ttot_Inv,tmp,nij,nPert,D,MP,lMax,MPq,CoC, &
                      EC,LP_context%ANr,Standard,nBas1,nTemp,LP_context%Q_Nuc,Bond_Threshold,Opt_Method,iPlot,iPrint,nSym)

! If the dear user has requested to get diffuse distributions
! associated to the multipoles, go here.
if (Diffuse(1)) then
  call mma_allocate(TP,nAtoms,label='ToPoint')
  call mma_allocate(MPp,nij,mElem,label='NotToPoint')
  MPp(:,:) = MP(:,:,1)
  call CoreToPoint(nAtoms,MPp,TP)
  LuYou = IsFreeUnit(81)
  call OpnFl('DIFFPR',LuYou,Exists)
  call Diff_MotherGoose(Diffuse,nAtoms,nBas1,MPp,nij,EC,LP_context%ANr,Ttot,Ttot_Inv,lMax,TP,dLimmo,Thrs1,Thrs2,nThrs,iPrint, &
                        ThrsMul,LuYou)
  close(LuYou)
  call mma_deallocate(TP)
  call mma_deallocate(MPp)
end if

call mma_deallocate(Origin)
call mma_deallocate(sq_mu)
call mma_deallocate(Ttot)
call mma_deallocate(Ttot_Inv)
call mma_deallocate(sq_temp)
call mma_deallocate(LP_context%center)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Pol,6,nij,label='Pol')
call mma_allocate(Cpl,6,nij,label='Cpl')
call mma_allocate(CplT,6,nij,label='CplT')

if (.not. NoField) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      D Y N A M I C   P R O P E R T I E S                           *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the fluctuating charges

  call Make_Fluctuating_Charges(nAtoms,LP_context%ANr,nij,nPert,MP,mElem,EC,Alpha)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Assemble the localized polarizabilities

  call Dynamic_Properties(tmp,nAtoms,MP,nij,nPert,mElem,Delta,EC,Pol,LP_context%ANr,Bond_Threshold,Cpl,CplT)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the properties

call mma_allocate(LbLCnt,nAtoms,label='LblCnt')
call mma_allocate(LbLCnt4,nAtoms,label='LblCnt4')
call Get_cArray('LP_L',LblCnt4,LenIn4*nAtoms)
do i=1,nAtoms
  LblCnt(i) = LblCnt4(i)(1:LenIn)
end do
call mma_deallocate(LblCnt4)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate arrays needed in Print_Local

nDim = nij*mElem
call mma_allocate(xMP,nDim,label='xMP')
call mma_allocate(xxMP,nDim,label='xxMP')
call mma_allocate(nxMP,nDim,label='nxMP')
call Print_Local(MP,nij,mElem,LP_context%C,nAtoms,CoC,LP_context%Q_Nuc,lMax,LblCnt,MPq,EC,Pol,NoField,tmp,xMP,xxMP,nxMP, &
                 LP_context%ANr,nOcOb,Energy_Without_FFPT,Ene_Occ,MpProp_Level,Bond_Threshold,CplT,LIonize)

call mma_deallocate(LblCnt)
call mma_deallocate(Ene_Occ)
call mma_deallocate(xMP)
call mma_deallocate(xxMP)
call mma_deallocate(nxMP)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Pol)
call mma_deallocate(Cpl)
call mma_deallocate(CplT)

call mma_deallocate(MPq)
call mma_deallocate(EC)
call mma_deallocate(MP)
call mma_deallocate(tmp)
call mma_deallocate(LP_context%Q_Nuc)
call mma_deallocate(LP_context%ANr)
call mma_deallocate(LP_context%C)
call mma_deallocate(LP_context%P)
call mma_deallocate(LP_context%PInv)

! Set flag on runfile to identify that it can be used to restart a
! LoProp calculation in the future.

if (.not. Restart) then
  LoProp_Mode = 1
  if (NoField) then
    LoProp_Mode = 2
  end if
  call Put_iScalar('LoProp Restart',LoProp_Mode)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Cleanup so that finish will not scream.

!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = 0
!                                                                      *
!***********************************************************************
!                                                                      *
end subroutine Polar
