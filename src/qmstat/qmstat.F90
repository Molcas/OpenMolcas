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
implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension Coord(MxAt*3)
dimension nBas(1), nBas_C(1), nCnC_C(MxBasC)
!******JoseMEP New variables to perform the MEP calculation
dimension Eint(MxQCen,10), Poli(MxQCen,10)
dimension SumElcPot(MxQCen,10)
dimension PertElcInt(MxBas*(MxBas+1)/2)
dimension nOcc(MxBas), natyp(MxAt)
character Labjhr*4
!***************

! The journey begins. Set non-zero return code.

ireturn = 99

! If parallel compilation, terminate, but gracefully.

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  write(6,*)
  write(6,*) ' QMStat does not run in parallel!'
  write(6,*)
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
  Go To 666
end if
if (Anal) then
  call Analyze_Q(iQ_Atoms)
  Go To 666
end if

! Read in orbitals, basis functions, integrals and bla bla bla. This
! is the centre for communicating with the rest of Molcas.

!******JoseMEP*** Qfread is called with more variables to the MEP calculation
call Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C,nOcc,natyp,nntyp)
!*******

! The turning point for the single-point calculation.

nCalls = 0
3333 continue

! If user request a set of single point calcualtions, then go in to
! a separate routine and do a 'reintrepretation' of the input.

if (SingPoint) call SingP(nCalls,iQ_Atoms,ipStoreCoo,nPart2)

! Decide which type of calculation we are to run. At this stage only
! QM equilibration and QM production is available. Should probably
! be so in the future, hence all-classical should be removed
! at some stage. Also, ordinary sampfile analysis is performed
! separately.

if (Smeq .or. Smprod) then
  write(6,*)
  write(6,*) 'All classical simulation is currently not available.'
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

if (SingPoint) then
  if (nCalls < nPart2) then
    write(6,*)
    write(6,*)
    write(6,*)
    write(6,*) '     *********** << NEW CALCULATION >> ***********'
    write(6,'(A,I3,A)') 'Single point-calculation nr.',nCalls+1,' follows below.'
    write(6,*)
    write(6,*)
    Go To 3333
  else
    call GetMem('Store','Free','Real',ipStoreCoo,nPart2*nCent*3)
    write(6,*)
    write(6,*)
    write(6,*) ' Single-point calculation ended.'
    write(6,*)
  end if
end if

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

! Exit

666 continue

ireturn = 0

return

end subroutine Qmstat
