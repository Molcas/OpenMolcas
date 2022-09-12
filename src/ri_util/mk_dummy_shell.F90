
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2008, Roland Lindh                                     *
!***********************************************************************

subroutine Mk_Dummy_Shell()
!***********************************************************************
!                                                                      *
!     Add the final DUMMY SHELL!                                       *
!                                                                      *
! 2008 R. Lindh, Dept. of Theor. Chem., Univ. of Lund, Sweden          *
!***********************************************************************

use Basis_Info, only: dbsc, iCnttp_Dummy, Max_Shells, nCnttp, Shells
use Center_Info, only: dc, n_dc
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iShll, mdc, nCnt, nCntrc, nPrim

!                                                                      *
!***********************************************************************
!                                                                      *
iShll = S%Mx_Shll-1
mdc = dbsc(nCnttp)%mdci+dbsc(nCnttp)%nCntr
nCnttp = nCnttp+1
if (nCnttp > Mxdbsc) then
  call WarningMessage(2,'Mk_Dummy_Shell: Increase Mxdbsc')
  call Abend()
end if
dbsc(nCnttp)%iVal = iShll+1
dbsc(nCnttp)%nVal = 1
dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal

dbsc(nCnttp)%Bsl = '.....RI_Dummy'
dbsc(nCnttp)%AtmNr = 1
dbsc(nCnttp)%Aux = .true.
dbsc(nCnttp)%Charge = Zero

nPrim = 1
nCntrc = 1

iShll = iShll+1
Shells(iShll)%Aux = .true.
call mma_allocate(Shells(iShll)%Exp,nPrim,Label='ExpDummy')
Shells(iShll)%nExp = nPrim
Shells(iShll)%nBasis = nCntrc
Shells(iShll)%nBasis_c = nCntrc
! Exponent
Shells(iShll)%Exp(1) = Zero
! Coefficients
call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
Shells(iShll)%Cff_c(1,1,1) = One
Shells(iShll)%Cff_c(1,1,2) = One
Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)

Shells(iShll)%Transf = .false.
Shells(iShll)%Prjct = .false.

! The coordinates

nCnt = 1
n_dc = max(mdc+nCnt,n_dc)
if (mdc+nCnt > MxAtom) then
  call WarningMessage(2,'Mk_Dummy_Shell: Increase MxAtom')
  call Abend()
end if
dbsc(nCnttp)%mdci = mdc
dc(mdc+nCnt)%LblCnt = 'Origin'
if (mdc+nCnt > 1) call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)
call mma_allocate(dbsc(nCnttp)%Coor_Hidden,3,1,Label='dbsc:C')
dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
dbsc(nCnttp)%Coor(1:3,1:1) = Zero
dbsc(nCnttp)%nCntr = nCnt
mdc = mdc+nCnt
!                                                                      *
!***********************************************************************
!                                                                      *
S%Mx_Shll = iShll+1
Max_Shells = S%Mx_Shll
S%Mx_mdc = mdc

if (iCnttp_Dummy /= 0) then
  write(u6,*) 'Mk_dummy_shell: iCnttp_Dummy'
  call Abend()
end if
iCnttp_Dummy = nCnttp
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_Dummy_Shell
