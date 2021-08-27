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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: ireturn
#include "itmax.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
integer(kind=iwp), parameter :: nElem=(iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6
integer(kind=iwp) :: i, ip_ANr, ip_Center, ip_D(0:6), ip_mu(0:nElem-1), ip_sq_mu(0:nElem-1), ip_Type, ipC, iPert, iPL, iPlot, ipP, &
                     ipPInv, ipQ_Nuc, iPrint, lMax, mElem, nAtoms, nBas(8), nBas1, nBas2, nBasMax, nij, nOrb(8), nPert, nSize, &
                     nStateF, nStateI, nSym, nTemp
real(kind=wp) :: Bond_Threshold, CoC(3), Origin(3,0:iTabMx), SubScale
logical(kind=iwp) :: lSave, NoField, PrintDen, Restart, Standard, SubtractDen, TDensity, UserDen, Utility
character(len=LenIn4) :: LblCnt(MxAtom)
character(len=12) :: Opt_Method
real(kind=wp), allocatable :: EC(:,:), MP(:,:,:), MPq(:), tmp(:), sq_temp(:), Ttot(:,:), Ttot_Inv(:,:)
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

call Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,ipC,ipQ_Nuc,ip_ANr,ip_Type,ip_Center,nSize,nBas1,nBas2,nBasMax,ipP,ipPInv)
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

call Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,iWork(ip_Center),iWork(ip_Type),nBas1,nBas2,nSym,nBasMax,ipPInv,Restart)

call Free_iWork(ip_type)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read in the multipole moment integrals once and for all.

lMax = 0   ! do only charges
mElem = (lMax*(lMax**2+6*lMax+11)+6)/6

nTemp = nBas1**2
call mma_allocate(tmp,nTemp,label='tmp')

call mma_allocate(MPq,mElem,label='MPq')
call Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Ttot,tmp,Origin,MPq,mElem,nBas1,nBas2,nBasMax,nTemp,nSym,ipPInv,Restart,Utility)
call mma_deallocate(Ttot)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the 1-particle density matrix

iPert = 0
call Get_Density_Matrix(ip_D(0),nBas1,nBas2,nBasMax,nBas,nSym,ipP,UserDen,PrintDen,SubtractDen,SubScale,Work(ipQ_Nuc),nAtoms, &
                        iPert,Restart,Utility,TDensity,nStateI,nStateF)

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
call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,sq_temp,Origin,iWork(ip_center),Ttot_Inv,tmp,nij,nPert,ip_D,MP,lMax,MPq,CoC, &
                      EC,iWork(ip_ANr),Standard,nBas1,nTemp,Work(ipQ_Nuc),Bond_Threshold,Opt_Method,iPlot,iPrint,nSym)

do i=mElem,1,-1
  call Free_Work(ip_sq_mu(i-1))
end do
call mma_deallocate(Ttot_Inv)
call mma_deallocate(sq_temp)
call Free_iWork(ip_center)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the properties

call Get_cArray('LP_L',LblCnt,LenIn4*nAtoms)
call LoProp_Print(MP,nij,nElem,nAtoms,Work(ipQ_Nuc),LblCnt,lSave)
if (iPL >= 2) then
  call CollapseOutput(0,'   Static properties:')
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Work(ipQ_Nuc)
call mma_deallocate(MPq)
call mma_deallocate(EC)
call mma_deallocate(MP)
call mma_deallocate(tmp)
call Free_iWork(ip_ANr)
call Free_Work(ipC)
if (nSym /= 1) then
  call Free_Work(ipP)
  call Free_Work(ipPInv)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = 0
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine LoProp
