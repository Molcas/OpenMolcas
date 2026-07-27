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
<<<<<<< HEAD
subroutine SG_setup_RASSI(nSym,nActEl,iSpin,SGS,CIS)

use Molcas, only: MxLev
use sguga, only: CIStruct, SGStruct, SG_Init
use rassi_aux, only: Level
use rassi_data, only: NASH
use definitions, only: iwp

integer(kind=iwp), intent(in):: nSym,nActEl,iSpin
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS

integer(kind=iwp) :: nLev, ISYM, IT, ILEV, ISM(MxLev), L2Act(MxLev), iq
=======

subroutine SG_setup_RASSI(nSym,nActEl,iSpin,SGS,CIS,EXS)

use Molcas, only: MxLev
use sguga, only: SGStruct, CIStruct, EXStruct, SG_Init
use rassi_aux, only: Level
use rassi_data, only: NASH
use rasdef, only: nRas, nRasEl, nRsPrt
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: ILEV, iq, ISM(MxLev), ISYM, IT, L2Act(MxLev), nLev
>>>>>>> upstream-openmolcas/master

nLev = 0
do ISYM=1,NSYM
  do IT=1,NASH(ISYM)
    nLev = nLev+1
<<<<<<< HEAD
    ILEV = LEVEL(nLev)
=======
    ILEV = Level(nLev)
>>>>>>> upstream-openmolcas/master
    ISM(ILEV) = ISYM
  end do
end do

<<<<<<< HEAD
L2Act(1:MxLev)=[(iq,iq=1,MxLev)]

Call SG_Init(nSym,nActEl,iSpin,SGS,CIS,                &
             xLevel=Level,xL2Act=L2Act,                &
             xNLEV=nLev,xNSM=ISM)

End subroutine SG_setup_RASSI
=======
L2Act(1:MxLev) = [(iq,iq=1,MxLev)]

call SG_Init(nSym,nActEl,iSpin,SGS,CIS,                    &
             nRas,nRasEl,nRsPrt,                           &
             EXS,                                          &
             xLevel=Level,xL2Act=L2Act,xNLEV=nLev,xNSM=ISM)

end subroutine SG_setup_RASSI
>>>>>>> upstream-openmolcas/master
