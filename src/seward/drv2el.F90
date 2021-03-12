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
! Copyright (C) 1990,1991,1993,1998, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El(Integral_WrOut,ThrAO)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!***********************************************************************

use iSD_data
use Basis_Info, only: dbsc
use Real_Info, only: CutInt

implicit real*8(A-H,O-Z)
external Integral_WrOut, Rsv_GTList
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
parameter(nTInt=1,mDens=1)
real*8 Dens(mDens), Fock(mDens), TInt(nTInt)
integer iTOffs(8,8,8)
logical Verbose, Indexation, FreeK2, W2Disc, PreSch, DoIntegrals, DoFock, DoGrad, FckNoClmb, FckNoExch, Rsv_GTList, Triangular
character*72 SLine
real*8, dimension(:,:), allocatable :: TMax
integer, dimension(:,:), allocatable :: Pair_Index
dimension ExFac(1), FckNoClmb(1), FckNoExch(1)

!                                                                      *
!***********************************************************************
!                                                                      *
SLine = 'Computing 2-electron integrals'
call StatusLine(' Seward:',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
ExFac = One
Nr_Dens = 1
DoIntegrals = .true.
DoFock = .false.
DoGrad = .false.
FckNoClmb = .false.
FckNoExch = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

Indexation = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
Thize = Zero               ! Not used for conventional integrals
PreSch = .true.            ! Not used for conventional integrals

Disc = Zero
Dix_Mx = Zero
TskHi = Zero
TskLw = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal)
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Pair_Index,2,nSkal*(nSkal+1)/2)
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = iS
      Pair_Index(2,nij) = jS
    end if
  end do
end do
P_Eff = dble(nij)
!                                                                      *
!***********************************************************************
!                                                                      *
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList
call Init_GTList
iOpt = 0

PP_Eff = P_Eff**2
PP_Eff_delta = 0.10d0*PP_Eff
PP_Count = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

! big loop over individual tasks distributed over individual nodes

10 continue
! make reservations of a tesk in global task list and get task range
! in return. Function will be false if no more tasks to execute.
if (.not. Rsv_GTlist(TskLw,TskHi,iOpt,W2Disc)) Go To 11
W2Disc = .false.

! Now do a quadruple loop over shells

ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
iS = Pair_Index(1,ijS)
jS = Pair_Index(2,ijS)
klS = int(TskLw-dble(ijS)*(dble(ijS)-One)/Two)
kS = Pair_Index(1,klS)
lS = Pair_Index(2,klS)
Count = TskLw

if (Count-TskHi > 1.0D-10) Go To 12
13 continue

! Logic to avoid computing integrals in a mixed muonic and
! electronic basis.

iCnttp = iSD(13,iS)
jCnttp = iSD(13,jS)
if (dbsc(iCnttp)%fMass /= dbsc(jCnttp)%fMass) Go To 14
kCnttp = iSD(13,kS)
lCnttp = iSD(13,lS)
if (dbsc(kCnttp)%fMass /= dbsc(lCnttp)%fMass) Go To 14

S_Eff = dble(ijS)
T_Eff = dble(klS)
ST_Eff = S_Eff*(S_Eff-One)/2d0+T_Eff
if (ST_Eff >= PP_Count) then
  write(SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',ST_Eff/PP_Eff*100d0,'% done so far.'
  call StatusLine(' Seward:',SLine)
  PP_Count = PP_Count+PP_Eff_delta
end if

Aint = TMax(iS,jS)*TMax(kS,lS)
if (AInt < CutInt) Go To 14
! from Dens are dummy arguments
call Eval_Ints_New_Inner(iS,jS,kS,lS,TInt,nTInt,iTOffs,Integral_WrOut,Dens,Fock,mDens,ExFac,Nr_Dens,FckNoClmb,FckNoExch,Thize, &
                         W2Disc,PreSch,Dix_Mx,Disc,Count,DoIntegrals,DoFock)
14 continue
Count = Count+One
if (Count-TskHi > 1.0D-10) Go To 12
klS = klS+1
if (klS > ijS) then
  ijS = ijS+1
  klS = 1
end if
iS = Pair_Index(1,ijS)
jS = Pair_Index(2,ijS)
kS = Pair_Index(1,klS)
lS = Pair_Index(2,klS)
Go To 13

! Task endpoint
12 continue

! Use a time slot to save the number of tasks and shell
! quadrupltes process by an individual node
call SavStat(1,One,'+')
call SavStat(2,TskHi-TskLw+One,'+')
Go To 10
11 continue
! End of big task loop
call CWTime(TCpu2,TWall2)
call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_GTList
call Free_PPList
call Free_TList

call mma_deallocate(Pair_Index)
call mma_deallocate(TMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
call Free_iSD()

return

end subroutine Drv2El
