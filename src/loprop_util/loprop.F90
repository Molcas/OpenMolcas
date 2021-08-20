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

implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
parameter(nElem=(iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6)
#include "WrkSpc.fh"
#include "real.fh"
real*8 Origin(3,0:iTabMx), CoC(3)
integer nBas(8), ip_mu(0:nElem-1), nOrb(8), ip_sq_mu(0:nElem-1), ip_D(0:6)
logical NoField, Standard, Utility, UserDen, PrintDen, SubtractDen
logical Restart, TDensity, lSave, Reduce_Prt
external Reduce_Prt
character*(LENIN4) LblCnt(MxAtom)
character*12 Opt_Method

!                                                                      *
!***********************************************************************
!                                                                      *
Utility = .true.
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
Bond_Threshold = -1.0d0
iPlot = 0
iPrint = 0
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

lMax = 0   ! do only charges
mElem = (lMax*(lMax**2+6*lMax+11)+6)/6

nTemp = nBas1**2
call Allocate_Work(ip_tmp,nTemp)

call Allocate_Work(ipMPq,mElem)
call Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Work(ip_Ttot),Work(ip_tmp),Origin,Work(ipMPq),mElem,nBas1,nBas2,nBasMax,nTemp, &
                        nSym,ipPInv,Restart,Utility)
call Free_Work(ip_Ttot)
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
nmu = nij*mElem*nPert
call Allocate_Work(ipMP,nmu)
call Allocate_Work(ip_sq_temp,nTemp)
call Allocate_Work(ip_EC,3*nij)

if (iPL >= 2) then
  write(6,*)
  call CollapseOutput(1,'   Static properties:')
  write(6,'(3X,A)') '   ------------------'
  write(6,*)
end if
call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,Work(ip_sq_temp),Origin,iWork(ip_center),Work(ip_Ttot_Inv),Work(ip_tmp),nij, &
                      nPert,ip_D,Work(ipMP),lMax,Work(ipMPq),CoC,Work(ip_EC),iWork(ip_ANr),Standard,nBas1,nTemp,Work(ipQ_Nuc), &
                      Bond_Threshold,Utility,Opt_Method,iPlot,iPrint,nSym)

do i=mElem,1,-1
  call Free_Work(ip_sq_mu(i-1))
end do
call Free_Work(ip_Ttot_Inv)
call Free_Work(ip_sq_temp)
call Free_iWork(ip_center)
!                                                                      *
!***********************************************************************
!                                                                      *
call Allocate_Work(ipPol,6*nij)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the properties

call Get_cArray('LP_L',LblCnt,(LENIN4)*nAtoms)
call LoProp_Print(Work(ipMP),nij,nElem,nAtoms,Work(ipQ_Nuc),LblCnt,lSave)
if (iPL >= 2) then
  call CollapseOutput(0,'   Static properties:')
  write(6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Work(ipPol)
call Free_Work(ipQ_Nuc)
call Free_Work(ipMPq)
call Free_Work(ip_EC)
call Free_Work(ipMP)
call Free_Work(ip_Tmp)
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
