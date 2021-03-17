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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "itmax.fh"
#include "Molcas.fh"
#include "timtra.fh"
#include "WrkSpc.fh"
integer(kind=iwp), parameter :: nElem = (iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6
integer(kind=iwp) :: i, i_f, ip_ANr, ip_Center, ip_D(0:6), ip_EC, ip_Ene_Occ, ip_h0, ip_mu(0:nElem-1), ip_sq_mu(0:nElem-1), &
                     ip_Ttot, ip_Ttot_Inv, ip_Type, ipC, iPert, iPlot, ipMPp, ipP, ipPInv, ipQ_Nuc, iPrint, ipXHole2, iTP, lMax, &
                     LoProp_Mode, LuYou, mElem, MpProp_Level, nAtoms, nBas(8), nBas1, nBas2, nBasMax, nDim, nij, nmu, nOcOb, &
                     nOrb(8), nPert, nSize, nStateF, nStateI, nSym, nTemp, nThrs
real(kind=wp) :: Alpha, Bond_Threshold, CoC(3), Delta, Dlt, dLimmo(2), dMolExpec, Energy_Without_FFPT, Ep, Origin(3,0:iTabMx), &
                 SubScale, Thrs1, Thrs2, ThrsMul
logical(kind=iwp) :: NoField, Standard, Utility, UserDen, PrintDen, SubtractDen, Restart, TDensity, XHole, Diffuse(3), Exists, &
                     LIonize
character(len=12) :: Opt_Method
character(len=LenIn), allocatable :: LblCnt(:)
character(len=LenIn4), allocatable :: LblCnt4(:)
real(kind=wp), allocatable :: Cpl(:,:), CplT(:,:), EC(:,:), MP(:), MPp(:), MPq(:), nxMP(:), Pol(:,:), sq_temp(:), tmp(:), TP(:), &
                              Ttot(:,:), Ttot_Inv(:,:), XHLoc2(:), xMP(:), xxMP(:)
integer(kind=iwp), external :: IsFreeUnit, ip_of_Work

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
Utility = .false.
ireturn = 99
!                                                                      *
!***********************************************************************
!                                                                      *
! Prelude

call Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,ipC,ipQ_Nuc,ip_ANr,ip_Type,ip_Center,nSize,nBas1,nBas2,nBasMax,ipP,ipPInv)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the input
!
! NoField is defaulted to true if symmetry is implied.

NoField = nSym /= 1
call Readin_polar(NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen,SubtractDen, &
                  SubScale,Restart,TDensity,nStateI,nStateF,XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,LIonize)
call InfoToMp(nSym,nBas,Energy_Without_FFPT,ip_Ene_Occ,nOcOb,UserDen,Restart)
!                                                                      *
!***********************************************************************
!                                                                      *
! Do the LoProp localization.

call mma_allocate(Ttot,nBas1,nBas1,label='Ttot')
call mma_allocate(Ttot_Inv,nBas1,nBas1,label='TtotInv')

call Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,iWork(ip_Center),iWork(ip_Type),nBas1,nBas2,nSym,nBasMax,ipPInv,Restart)

call Free_iWork(ip_type)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read in the multipole moment integrals once and for all.

call Get_iScalar('Highest Mltpl',lMax)
write(u6,'(A,I2)') ' Multipole moments will be processed up to order ',lMax
write(u6,*)
mElem = (lMax*(lMax**2+6*lMax+11)+6)/6

nTemp = nBas1**2
call mma_allocate(tmp,nTemp,label='tmp')

call mma_allocate(MPq,mElem,label='MPq')
call Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Ttot,tmp,Origin,MPq,mElem,nBas1,nBas2,nBasMax,nTemp,nSym,ipPInv,Restart,Utility)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the 1-particle density matrix

Dlt = -Delta
iPert = 0
call Get_Density_Matrix(ip_D(0),nBas1,nBas2,nBasMax,nBas,nSym,ipP,UserDen,PrintDen,SubtractDen,SubScale,Work(ipQ_Nuc),nAtoms, &
                        iPert,Restart,Utility,TDensity,nStateI,nStateF)

! If computing local xhole-dipole moments. Should come after
! get_density_matrix so modified densities are used.
if (XHole) call Compute_XHole_Int(nBas,nSym,ipXHole2,dMolExpec,nSize)

if (.not. NoField) then
  ! Read the one-electron hamiltonian.
  call Read_h0(nSize,nBas(1),ip_h0,Restart)
  do iPert=1,6
    i_f = (iPert+1)/2
    Dlt = -Dlt
    if ((.not. Restart) .and. (.not. UserDen)) then
      call Comp_F(Work(ip_h0),Work(ip_mu(i_f)),nBas(1),Dlt,Ep,Work(ip_mu(0)),CoC(i_f),Origin(i_f,1))
    end if
    call Get_Density_Matrix(ip_D(iPert),nBas1,nBas2,nBasMax,nBas,nSym,ipP,UserDen,PrintDen,SubtractDen,SubScale,Work(ipQ_Nuc), &
                            nAtoms,iPert,Restart,Utility,TDensity,nStateI,nStateF)
  end do
  call Free_Work(ip_h0)
end if
do i=mElem,1,-1
  call Free_Work(ip_mu(i-1))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!     S T A T I C   P R O P E R T I E S                                *
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute localized multipole moments

nPert = 2*3+1
if (NoField) nPert = 1
nij = (nAtoms*(nAtoms+1)/2)
nmu = nij*mElem*nPert
call mma_allocate(MP,nmu,label='MP')
call mma_allocate(sq_temp,nTemp,label='sq_temp')
call mma_allocate(EC,3,nij,label='EC')

call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,sq_temp,Origin,iWork(ip_center),Ttot_Inv,tmp,nij,nPert,ip_D,MP,lMax,MPq,CoC, &
                      EC,iWork(ip_ANr),Standard,nBas1,nTemp,Work(ipQ_Nuc),Bond_Threshold,Utility,Opt_Method,iPlot,iPrint,nSym)

!-- If XHole integrals are available, localize them. Most unfortunate,
!   the local_properties routine is focused on multipole moments,
!   hence we rather write a new routine for Xhole, than significantly
!   edit the local_properties routine.
if (XHole) then
  call mma_allocate(XHLoc2,nij,label='XHLoc2')
  call Local_Xhole(ipXHole2,dMolExpec,nAtoms,nBas1,nTemp,iWork(ip_center),Ttot,Ttot_Inv,Work(ipC),nij,EC,iWork(ip_ANr), &
                   Bond_Threshold,iPrint,XHLoc2)
else
  call mma_allocate(XHLoc2,0,label='XHLoc2')
end if

!-- If the dear user has requested to get diffuse distributions
!   associated to the multipoles, go here.
if (Diffuse(1)) then
  call mma_allocate(TP,nAtoms,label='ToPoint')
  call mma_allocate(MPp,nmu,label='NotToPoint')
  call dcopy_(nmu,MP,1,MPp,1)
  call CoreToPoint(nAtoms,MPp,TP)
  LuYou = IsFreeUnit(81)
  call OpnFl('DIFFPR',LuYou,Exists)
  iTP = ip_of_Work(TP(1))
  ipMPp = ip_of_Work(MPp(1))
  ip_Ttot = ip_of_Work(Ttot(1,1))
  ip_Ttot_Inv = ip_of_Work(Ttot_Inv(1,1))
  ip_EC = ip_of_Work(EC(1,1))
  call Diff_MotherGoose(Diffuse,nAtoms,nBas1,ipMPp,ipC,nij,ip_EC,ip_ANr,ip_Ttot,ip_Ttot_Inv,lMax,iTP,dLimmo,Thrs1,Thrs2,nThrs, &
                        iPrint,ThrsMul,LuYou)
  close(LuYou)
  call mma_deallocate(TP)
  call mma_deallocate(MPp)
end if

do i=mElem,1,-1
  call Free_Work(ip_sq_mu(i-1))
end do
call mma_deallocate(Ttot)
call mma_deallocate(Ttot_Inv)
call mma_deallocate(sq_temp)
call Free_iWork(ip_center)
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

  call Make_Fluctuating_Charges(nAtoms,iWork(ip_ANr),nij,nPert,MP,mElem,EC,Alpha)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Assemble the localized polarizabilities

  call Dynamic_Properties(tmp,nAtoms,MP,nij,nPert,mElem,Delta,EC,Pol,iWork(ip_ANr),Bond_Threshold,Cpl,CplT)

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
call Print_Local(MP,nij,mElem,Work(ipC),nAtoms,CoC,Work(ipQ_Nuc),lMax,LblCnt,MPq,EC,Pol,NoField,tmp,xMP,xxMP,nxMP,iWork(ip_ANr), &
                 nOcOb,Energy_Without_FFPT,Work(ip_Ene_Occ),MpProp_Level,Bond_Threshold,XHole,XHLoc2,dMolExpec,CplT,LIonize)

call mma_deallocate(LblCnt)
call Free_Work(ip_Ene_Occ)
call mma_deallocate(xMP)
call mma_deallocate(xxMP)
call mma_deallocate(nxMP)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Pol)
call mma_deallocate(Cpl)
call mma_deallocate(CplT)

call Free_Work(ipQ_Nuc)
call mma_deallocate(MPq)
call mma_deallocate(EC)
call mma_deallocate(MP)
call mma_deallocate(tmp)
call Free_iWork(ip_ANr)
call Free_Work(ipC)
call mma_deallocate(XHLoc2)
if (nSym /= 1) then
  call Free_Work(ipP)
  call Free_Work(ipPInv)
end if

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

nfld_tim = 0
nfld_stat = 0
!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = 0
!                                                                      *
!***********************************************************************
!                                                                      *
end subroutine Polar
