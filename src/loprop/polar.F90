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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "itmax.fh"
#include "Molcas.fh"
#include "timtra.fh"
#include "WrkSpc.fh"
integer(kind=iwp), parameter :: nElem = (iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6
integer(kind=iwp) :: i, i_f, ip_ANr, ip_Center, ip_D(0:6), ip_EC, ip_Ene_Occ, ip_h0, ip_mu(0:nElem-1), ip_sq_mu(0:nElem-1), &
                     ip_sq_temp, ip_tmp, ip_Ttot, ip_Ttot_Inv, ip_Type, ipC, ipCpl, ipCplT, iPert, iPlot, ipMP, ipMPp, ipMPq, &
                     ipnxMP, ipP, ipPInv, ipPol, ipQ_Nuc, iPrint, ipXHLoc2, ipXHole2, ipxMP, ipxxMP, iTP, lMax, LoProp_Mode, &
                     LuYou, mElem, MpProp_Level, nAtoms, nBas(8), nBas1, nBas2, nBasMax, nDim, nij, nmu, nOcOb, nOrb(8), nPert, &
                     nSize, nStateF, nStateI, nSym, nTemp, nThrs
real(kind=wp) :: Alpha, Bond_Threshold, CoC(3), Delta, Dlt, dLimmo(2), dMolExpec, Energy_Without_FFPT, Ep, Origin(3,0:iTabMx), &
                 SubScale, Thrs1, Thrs2, ThrsMul
logical(kind=iwp) :: NoField, Standard, Utility, UserDen, PrintDen, SubtractDen, Restart, TDensity, XHole, Diffuse(3), Exists, &
                     LIonize
character(len=LenIn) :: LblCnt(MxAtom)
character(len=LenIn4) :: LblCnt4(MxAtom)
character(len=12) :: Opt_Method
integer(kind=iwp), external :: IsFreeUnit

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

call GetMem('Ttot','Allo','Real',ip_Ttot,nBas1**2)
call GetMem('TtotInv','Allo','Real',ip_Ttot_Inv,nBas1**2)

call Localize_LoProp_Drv(Work(ip_Ttot),Work(ip_Ttot_Inv),nBas,iWork(ip_Center),iWork(ip_Type),nBas1,nBas2,nSym,nBasMax,ipPInv, &
                         Restart)

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
call Allocate_Work(ip_tmp,nTemp)

call Allocate_Work(ipMPq,mElem)
call Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Work(ip_Ttot),Work(ip_tmp),Origin,Work(ipMPq),mElem,nBas1,nBas2,nBasMax,nTemp, &
                        nSym,ipPInv,Restart,Utility)
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
call Allocate_Work(ipMP,nmu)
call Allocate_Work(ip_sq_temp,nTemp)
call Allocate_Work(ip_EC,3*nij)

call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,Work(ip_sq_temp),Origin,iWork(ip_center),Work(ip_Ttot_Inv),Work(ip_tmp),nij, &
                      nPert,ip_D,Work(ipMP),lMax,Work(ipMPq),CoC,Work(ip_EC),iWork(ip_ANr),Standard,nBas1,nTemp,Work(ipQ_Nuc), &
                      Bond_Threshold,Utility,Opt_Method,iPlot,iPrint,nSym)

!-- If XHole integrals are available, localize them. Most unfortunate,
!   the local_properties routine is focused on multipole moments,
!   hence we rather write a new routine for Xhole, than significantly
!   edit the local_properties routine.
if (XHole) then
  call Allocate_Work(ipXHLoc2,nij)
  call Local_Xhole(ipXHole2,dMolExpec,nAtoms,nBas1,nTemp,iWork(ip_center),Work(ip_Ttot),Work(ip_Ttot_Inv),Work(ipC),nij, &
                   Work(ip_EC),iWork(ip_ANr),Bond_Threshold,iPrint,ipXHLoc2)
  call Free_Work(ipXHole2)
else
  ipXHLoc2 = ip_Dummy
end if

!-- If the dear user has requested to get diffuse distributions
!   associated to the multipoles, go here.
if (Diffuse(1)) then
  call GetMem('ToPoint','Allo','Real',iTP,nAtoms)
  call GetMem('NotToPoint','Allo','Real',ipMPp,nmu)
  call dcopy_(nmu,Work(ipMP),1,Work(ipMPp),1)
  call CoreToPoint(nAtoms,ipMPp,iTP)
  LuYou = IsFreeUnit(81)
  call OpnFl('DIFFPR',LuYou,Exists)
  call Diff_MotherGoose(Diffuse,nAtoms,nBas1,ipMPp,ipC,nij,ip_EC,ip_ANr,ip_Ttot,ip_Ttot_Inv,lMax,iTP,dLimmo,Thrs1,Thrs2,nThrs, &
                        iPrint,ThrsMul,LuYou)
  close(LuYou)
  call GetMem('ToPoint','Free','Real',iTP,nAtoms)
  call GetMem('NotToPoint','Free','Real',ipMPp,nmu)
end if

do i=mElem,1,-1
  call Free_Work(ip_sq_mu(i-1))
end do
call Free_Work(ip_Ttot)
call Free_Work(ip_Ttot_Inv)
call Free_Work(ip_sq_temp)
call Free_iWork(ip_center)
!                                                                      *
!***********************************************************************
!                                                                      *
call Allocate_Work(ipPol,6*nij)
call Allocate_Work(ipCpl,6*nij)
call Allocate_Work(ipCplT,6*nij)

if (.not. NoField) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      D Y N A M I C   P R O P E R T I E S                           *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the fluctuating charges

  call Make_Fluctuating_Charges(nAtoms,iWork(ip_ANr),nij,nPert,Work(ipMP),mElem,Work(ip_EC),Alpha)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Assemble the localized polarizabilities

  call Dynamic_Properties(Work(ip_tmp),nAtoms,Work(ipMP),nij,nPert,mElem,Delta,Work(ip_EC),Work(ipPol),iWork(ip_ANr), &
                          Bond_Threshold,Work(ipCpl),Work(ipCplT))


end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the properties

call Get_cArray('LP_L',LblCnt4,(LENIN4)*nAtoms)
do i=1,nAtoms
  LblCnt(i)(1:LENIN) = LblCnt4(i)(1:LENIN)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate arrays needed in Print_Local

nDim = nij*mElem
call Allocate_Work(ipxMP,nDim)
call Allocate_Work(ipxxMP,nDim)
call Allocate_Work(ipnxMP,nDim)
call Print_Local(Work(ipMP),nij,mElem,Work(ipC),nAtoms,CoC,Work(ipQ_Nuc),lMax,LblCnt,Work(ipMPq),Work(ip_EC),Work(ipPol),NoField, &
                 Work(ip_Tmp),Work(ipxMP),Work(ipxxMP),Work(ipnxMP),iWork(ip_ANr),nOcOb,Energy_Without_FFPT,ip_Ene_Occ, &
                 MpProp_Level,Bond_Threshold,XHole,Work(ipXHLoc2),dMolExpec,Work(ipCpl),Work(ipCplT),LIonize)

call Free_Work(ipxMP)
call Free_Work(ipxxMP)
call Free_Work(ipnxMP)
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Work(ipPol)
call Free_Work(ipCpl)
call Free_Work(ipCplT)

call Free_Work(ipQ_Nuc)
call Free_Work(ipMPq)
call Free_Work(ip_EC)
call Free_Work(ipMP)
call Free_Work(ip_Tmp)
call Free_iWork(ip_ANr)
call Free_Work(ipC)
if (Xhole) call Free_Work(ipXHLoc2)
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
