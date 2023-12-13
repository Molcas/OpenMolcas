!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Qmstat(ireturn)

use qmstat_global, only: Anal, AvRed, BasOri, BigT, CasOri, ChaNuc, Cordst, DipIm, EdSt, info_atom, iQang, iQn, iWoGehenC, &
                         iWoGehenQ, lExtr, MoAveRed, mPrimus, MxPut, nBA_C, nBA_Q, nBonA_C, nBonA_Q, nCBoA_C, nCBoA_Q, nCent, &
                         nCnC_C, nPart, nPol, nPrimus, OldGeo, PertNElcInt, QIm, QImp, Qmeq, QmProd, Qmstat_end, QmType, SavOri, &
                         SingPoint, Sqrs, SupM, Trans, Udisp, V1
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: iCi, iQ_Atoms, nAtomsCC, nBas(1), nBas_C(1), nCalls, NCountField, nntyp
character(len=4) :: Labjhr
integer(kind=iwp), allocatable :: natyp(:), nOcc(:)
real(kind=wp), allocatable :: Coord(:,:), Eint(:,:), PertElcInt(:), Poli(:,:), StoreCoo(:,:,:), SumElcPot(:,:)
#include "warnings.h"
!******JoseMEP New variables to perform the MEP calculation
!Eint, Poli, SumElcPot, PertElcInt, nOcc, natyp, Labjhr

! The journey begins. Set non-zero return code.

ireturn = 99

! If parallel compilation, terminate, but gracefully.

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  write(u6,*)
  write(u6,*) ' QMStat does not run in parallel!'
  write(u6,*)
  call Quit(_RC_NOT_AVAILABLE_)
end if
#endif

! Set defaults, zeros and initial values.

call Qmstat_init()

! Make discrete banner.

! Read ('process') input. To do that we need the number of atoms in
! the QM-region, therefore that initial call to the RUNFILE.

call Get_iScalar('Unique atoms',iQ_Atoms)
call mma_allocate(info_atom,iQ_Atoms,label='info_atom')
call mma_allocate(ChaNuc,iQ_Atoms,label='ChaNuc')
call mma_allocate(Udisp,2,iQ_Atoms,label='Udisp')
uDisp(:,:) = Zero
call Get_Qmstat_Input(iQ_Atoms)

! If only the startfile is to be edited or the sampfile analyzed, go here, then terminate.

if (EdSt) then

  call EditStart()

else if (Anal) then

  call Analyze_Q(iQ_Atoms)

else

  ! Read in orbitals, basis functions, integrals and bla bla bla. This
  ! is the centre for communicating with the rest of Molcas.

  call mma_allocate(Coord,3,iQ_Atoms,label='Coord')
  call mma_allocate(nOcc,iQ_Atoms,label='nOcc')
  call mma_allocate(natyp,iQ_Atoms,label='natyp')

  !******JoseMEP*** Qfread is called with more variables to the MEP calculation
  call Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nOcc,natyp,nntyp)
  !*******

  call mma_allocate(PertElcInt,nTri_Elem(nBas(1)),label='PertElcInt')
  call mma_allocate(Eint,nTri_Elem(iQ_Atoms),10,label='Eint')
  call mma_allocate(Poli,nTri_Elem(iQ_Atoms),10,label='Poli')
  call mma_allocate(SumElcPot,nTri_Elem(iQ_Atoms),10,label='SumElcPot')
  call mma_allocate(OldGeo,3,size(Cordst,2),label='OldGeo')
  call mma_allocate(DipIm,3,MxPut*max(nCent,nPol),label='DipIm')
  call mma_allocate(QIm,MxPut*nCent,label='QIm')
  call mma_allocate(QImp,MxPut*max(nCent,nPol),label='QImp')
  call mma_allocate(Sqrs,MxPut*nCent,label='Sqrs')
  call mma_allocate(PertNElcInt,nTri_Elem(nBas(1)),label='PertNElcInt')
  call mma_allocate(CasOri,3,size(SavOri,2),label='CasOri')

  ! The turning point for the single-point calculation.

  nCalls = 0
  do

    ! If user request a set of single point calculations, then go in to
    ! a separate routine and do a 'reintrepretation' of the input.

    if (SingPoint) then
      if (nCalls == 0) call mma_allocate(StoreCoo,3,nCent,nPart,label='Store')
      call SingP(nCalls,iQ_Atoms,StoreCoo)
    end if

    ! Decide which type of calculation we are to run. At this stage only
    ! QM equilibration and QM production is available.
    ! Also, ordinary sampfile analysis is performed separately.

    if (Qmeq .or. Qmprod) then  !Qmeq=.true. is default option.
      if (QmType(1:3) == 'SCF') then
        call EqScf(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C)
      else if (QmType(1:4) == 'RASS') then
        call EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C)
      end if
    end if

    !********JoseMEP
    ! If MEP option true. Calculate the Interaction energy coming from the
    ! Interaction with Mean Elec. Potential, Field and Field-Gradient
    ! E=q*Pot-dipol*field+quad*field-grad
    ! Also is included the non_Electr perturbation
    ! For a quantum AO i nonE=d*Sum_k <Xi|Chi_k>*Sum_k <Xi|Chi_k>*E_k
    ! These Energies are added to the SEWARD One-electron file
    ! as a perturbation to the One-e Hamiltoniam

    if (lExtr(8)) then
      Labjhr = 'Pert'
      call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot,NCountField,PertElcInt,iQ_Atoms,nBas(1),nOcc,natyp,nntyp)
    end if
    !********

    ! If we are doing single point calculations, then jump back up.

    if (.not. SingPoint) exit
    if (nCalls < size(StoreCoo,3)) then
      write(u6,*)
      write(u6,*)
      write(u6,*)
      write(u6,*) '     *********** << NEW CALCULATION >> ***********'
      write(u6,'(A,I3,A)') 'Single point-calculation nr.',nCalls+1,' follows below.'
      write(u6,*)
      write(u6,*)
    else
      call mma_deallocate(StoreCoo)
      write(u6,*)
      write(u6,*)
      write(u6,*) ' Single-point calculation ended.'
      write(u6,*)
      exit
    end if
  end do

  call mma_deallocate(Coord)
  call mma_deallocate(nOcc)
  call mma_deallocate(natyp)
  call mma_deallocate(Eint)
  call mma_deallocate(Poli)
  call mma_deallocate(SumElcPot)
  call mma_deallocate(PertElcInt)
  call mma_deallocate(OldGeo)
  call mma_deallocate(DipIm)
  call mma_deallocate(QIm)
  call mma_deallocate(QImp)
  call mma_deallocate(Sqrs)
  call mma_deallocate(PertNElcInt)
  call mma_deallocate(CasOri)
  call mma_deallocate(nBA_C)
  call mma_deallocate(nBA_Q)
  call mma_deallocate(nBonA_C)
  call mma_deallocate(nBonA_Q)
  call mma_deallocate(nCBoA_C)
  call mma_deallocate(nCBoA_Q)
  call mma_deallocate(nCnC_C)
  call mma_deallocate(iWoGehenC)
  call mma_deallocate(iWoGehenQ)
  call mma_deallocate(iQn)
  call mma_deallocate(iQang)
  call mma_deallocate(BasOri)
  call mma_deallocate(SavOri)
  call mma_deallocate(mPrimus)
  call mma_deallocate(nPrimus)
  call mma_deallocate(Trans)

  ! Deallocations of stuff from qfread to simulations. These are unique for the QM-method.

  if (QmType(1:4) == 'RASS') then
    call mma_deallocate(BigT)
    if (MoAveRed) call mma_deallocate(AvRed)
  else if (QmType(1:3) == 'SCF') then
    call mma_deallocate(SupM)
    call mma_deallocate(V1)
  end if

end if

call mma_deallocate(info_atom)
call mma_deallocate(ChaNuc)
call mma_deallocate(Udisp)

! Exit

call Qmstat_end()

ireturn = 0

return

end subroutine Qmstat
