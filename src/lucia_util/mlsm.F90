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

subroutine MLSM(IML,IPARI,ISM,TYP,IWAY)
! Transfer between ML,IPARI notation and compound notation ISM
!
! IWAY = 1 : IML,IPARI => Compound
! IWAY = 2 : IML,IPARI <= Compound
!
! TYP : 'SX','OB','ST','DX','CI'

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IML, IPARI, ISM, IWAY
character(len=2) :: TYP
integer(kind=iwp) :: MNML, NML_, NTEST
integer(kind=iwp), parameter :: MNMLCI = 0, MNMLDX = 0, MNMLOB = 0, MNMLST = 0, MNMLSX = 0, NMLCI = 0, NMLDX = 0, NMLOB = 0, &
                                NMLST = 0, NMLSX = 0

select case (TYP)
  case ('OB')
    NML_ = NMLOB
    MNML = MNMLOB
  case ('SX')
    NML_ = NMLSX
    MNML = MNMLSX
  case ('DX')
    NML_ = NMLDX
    MNML = MNMLDX
  case ('ST')
    NML_ = NMLST
    MNML = MNMLST
  case ('CI')
    NML_ = NMLCI
    MNML = MNMLCI
  case default
    call Abend()
end select

if (IWAY == 1) then
  !ISM = (IPARI-1)*NML_+MNML-1
  ISM = (IPARI-1)*NML_+IML-MNML+1
else if (IWAY == 2) then
  if (ISM > NML_) then
    IPARI = 2
    IML = ISM-NML_+MNML-1
  else
    IPARI = 1
    IML = ISM+MNML-1
  end if
else
  write(u6,*) ' Error in MLSM, IWAY = ',IWAY
  write(u6,*) ' MLSM stop !!!'
  !stop 20
  call SYSABENDMSG('lucia_util/mlsm','Internal error','')
end if

NTEST = 0
if (NTEST /= 0) then
  write(u6,'(A,A)') ' MLSM speaking, type= ',TYP
  write(u6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
end if

end subroutine MLSM
