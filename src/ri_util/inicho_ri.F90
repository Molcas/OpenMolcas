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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  IniCho_RI
!
!> @brief
!>   Initialize Cholesky environment for RI calculations
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Initialize Cholesky environment for RI calculations.
!>
!> @note
!> Needs a call to ::SetUp_Ints with indexation turned on
!>
!> @param[in] nSkal    The number of shells (excl. aux. basis)
!> @param[in] nVec_Aux Number of aux. basis vectors per irrep
!> @param[in] nIrrep   Number of irreps
!> @param[in] iTOffs   Offset vector
!> @param[in] iShij    Index vector of shell pairs
!> @param[in] nShij    Number of shell pairs
!***********************************************************************

subroutine IniCho_RI(nSkal,nVec_Aux,nIrrep,iTOffs,iShij,nShij)

use Index_Functions, only: iTri
use RICD_Info, only: Thrshld_CD
use Para_Info, only: Is_Real_Par
use Cholesky, only: Cho_DecAlg, CHO_FAKE_PAR, Cho_Real_Par, IfcSew, InfRed, InfVec, IPRINT, iSP2F, MaxRed, MaxVec, nnShl, NumCho, &
                    nShell, nSym, RUN_EXTERNAL, RUN_MODE, ThrCom, XnPass
#ifdef _MOLCAS_MPP_
use Cholesky, only: myNumCho, NumCho_G
#endif
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSkal, nIrrep, nVec_Aux(0:nIrrep-1), iTOffs(3,nIrrep), nShij, iShij(2,nShij)
#ifdef _MOLCAS_MPP_
#endif
integer(kind=iwp) :: iDummy, ijS, iSym, iVec, LuOut
logical(kind=iwp) :: Alloc_Bkm, SetDefaultsOnly, Skip_PreScreen

! Set defaults for those parameters that can normally be changed
! through user input to the Cholesky decomposition.
! --------------------------------------------------------------

SetDefaultsOnly = .true.
iDummy = -1
LuOut = u6
call Cho_Inp(SetDefaultsOnly,iDummy,LuOut)

! Reset Cholesky Threshold for RI
! -------------------------------
ThrCom = Thrshld_CD

! Reset parallel config.
! ----------------------

CHO_FAKE_PAR = .false.
Cho_Real_Par = Is_Real_Par() .and. (.not. CHO_FAKE_PAR)

! Set run mode to "external" (should be irrelevant for RI).
! ---------------------------------------------------------

RUN_MODE = RUN_EXTERNAL

! Silence the Cholesky routines.
! ------------------------------

iPrint = 0

! Set number of shells (excl. aux. basis) in Cholesky
! ---------------------------------------------------

nShell = nSkal

! To avoid unnecessary allocations of shell-pair-to-reduced-set
! maps, set decomposition algorithm to 1 ("one-step") and the Seward
! interface to "1" (full shell quadruple storage). Both values are,
! obviously, irrelevant for RI. In parallel runs, use default values
! to avoid warnings being printed in Cho_P_Check.
! ------------------------------------------------------------------

if (Is_Real_Par()) then
  Cho_DecAlg = 4
  IfcSew = 2
else
  Cho_DecAlg = 1
  IfcSew = 1
end if

! Change MaxRed to 1 (all vectors have identical dimension, namely
! full => only 1 reduced set).
! ----------------------------------------------------------------

MaxRed = 1

! Set MaxVec to the largest number of vectors (= number of linearly
! independent auxiliary basis functions). In this way we avoid
! allocating more memory for InfVec than needed.
! -----------------------------------------------------------------

MaxVec = nVec_Aux(0)
do iSym=1,nIrrep-1
  MaxVec = max(MaxVec,nVec_Aux(iSym))
end do

! Other initializations. Most importantly, allocate InfRed and
! InfVec arrays (defined in Cholesky module).
! We skip diagonal prescreening, as it has already been done.
! Instead, allocate and set the mapping from reduced to full shell
! pairs here.
! ----------------------------------------------------------------

nnShl = nShij
call mma_allocate(iSP2F,nnShl,Label='iSP2F')
do ijS=1,nnShl
  iSP2F(ijS) = iTri(iShij(1,ijS),iShij(2,ijS))
end do
Skip_PreScreen = .true.
Alloc_Bkm = .false.
call Cho_Init(Skip_PreScreen,Alloc_Bkm)

! Set number of vectors equal to the number of lin. indep. auxiliary
! basis functions.
! ------------------------------------------------------------------

NumCho(1:nSym) = nVec_Aux(0:nSym-1)
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  NumCho_g(1:nSym) = NumCho(1:nSym)
  myNumCho(1:nSym) = 0
end if
#endif

! Do allocations that are normally done during or after the
! computation of the diagonal (since the dimension of the 1st
! reduced set is unknown until the screened diagonal is known).
! -------------------------------------------------------------

call IniCho_RI_Xtras(iTOffs,nIrrep,iShij,nShij)

! Set start disk addresses.
! -------------------------

XnPass = 0 ! it should be zeroed in Cho_Inp, but just in case.
call Cho_SetAddr(InfRed,InfVec,MaxRed,MaxVec,size(InfVec,2),nSym)

! Set vector info.
! Parent diagonal is set equal to the vector number, parent pass
! (i.e. reduced set) to 1.
! --------------------------------------------------------------

do iSym=1,nSym
  do iVec=1,NumCho(iSym)
    call Cho_SetVecInf(iVec,iSym,iVec,1,1)
  end do
end do

return

end subroutine IniCho_RI
