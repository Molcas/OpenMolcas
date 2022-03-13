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

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ireturn
#include "maxi.fh"
#include "qminp.fh"
#include "warnings.h"
integer(kind=iwp) :: iBigDeAll, iCi, ip, ipStoreCoo, iQ_Atoms, iSupDeAll, iV1DeAll, nAtomsCC, natyp(MxAt), nBas(1), nBas_C(1), & !IFG
                     nCalls, nCnC_C(MxBasC), NCountField, nntyp, nOcc(MxBas), nPart2, nRM, nS, nSizeBig !IFG
real(kind=wp) :: Coord(MxAt*3), Eint(MxQCen,10), PertElcInt(nTri_Elem(MxBas)), Poli(MxQCen,10), SumElcPot(MxQCen,10) !IFG
character(len=4) :: Labjhr
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
call Get_Qmstat_Input(iQ_Atoms)

! If only the startfile is to be edited or the sampfile analyzed, go here, then terminate.

if (EdSt) then

  call EditStart()

else if (Anal) then

  call Analyze_Q(iQ_Atoms)

else

  ! Read in orbitals, basis functions, integrals and bla bla bla. This
  ! is the centre for communicating with the rest of Molcas.

  !******JoseMEP*** Qfread is called with more variables to the MEP calculation
  call Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C,nOcc,natyp,nntyp)
  !*******

  ! The turning point for the single-point calculation.

  nCalls = 0
  do

    ! If user request a set of single point calcualtions, then go in to
    ! a separate routine and do a 'reintrepretation' of the input.

    if (SingPoint) call SingP(nCalls,iQ_Atoms,ipStoreCoo,nPart2)

    ! Decide which type of calculation we are to run. At this stage only
    ! QM equilibration and QM production is available. Should probably
    ! be so in the future, hence all-classical should be removed
    ! at some stage. Also, ordinary sampfile analysis is performed
    ! separately.

    if (Smeq .or. Smprod) then
      write(u6,*)
      write(u6,*) 'All classical simulation is currently not available.'
      call Quit(_RC_GENERAL_ERROR_)
    else if (Qmeq .or. Qmprod) then  !Qmeq=.true. is default option.
      if (QmType(1:3) == 'SCF') then
        call EqScf(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C,iSupDeAll,iV1DeAll)
      else if (QmType(1:4) == 'RASS') then
        call EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C,iBigDeAll,nSizeBig,ip,nRM)
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
      call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot,NCountField,PertElcInt,iQ_Atoms,nBas(1),nOcc(1),natyp(1),nntyp)
    end if
    !********

    ! If we are doing single point calculations, then jump back up.

    if (.not. SingPoint) exit
    if (nCalls < nPart2) then
      write(u6,*)
      write(u6,*)
      write(u6,*)
      write(u6,*) '     *********** << NEW CALCULATION >> ***********'
      write(u6,'(A,I3,A)') 'Single point-calculation nr.',nCalls+1,' follows below.'
      write(u6,*)
      write(u6,*)
    else
      call GetMem('Store','Free','Real',ipStoreCoo,nPart2*nCent*3)
      write(u6,*)
      write(u6,*)
      write(u6,*) ' Single-point calculation ended.'
      write(u6,*)
      exit
    end if
  end do

  ! Deallocations of stuff from qfread to simulations. These are unique for the QM-method.

  if (QmType(1:4) == 'RASS') then
    call GetMem('ALLES','Free','Real',iBigDeAll,nSizeBig)
    if (MoAveRed) then
      call GetMem('UncleMoe','Free','Real',ip,nBas(1)*nRM)
    end if
  else if (QmType(1:3) == 'SCF') then
    nS = iOrb(1)*(iOrb(1)+1)/2
    call GetMem('SUPER','Free','Real',iSupDeAll,nS**2)
    call GetMem('OrbCoeffQ','Free','Real',iV1DeAll,iOrb(1)*nBas(1))
  end if

end if

! Exit

ireturn = 0

return

end subroutine Qmstat
