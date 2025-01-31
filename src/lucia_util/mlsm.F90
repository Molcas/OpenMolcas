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

subroutine MLSM(IML,IPARI,ISM,type,IWAY)
! Transfer between ML,IPARI notation and compound notation ISM
!
! IWAY = 1 : IML,IPARI => Compound
! IWAY = 2 : IML,IPARI <= Compound
!
! TYPE : 'SX','OB','ST','DX','CI'

implicit none
integer IML, IPARI, ISM, IWAY
character(len=2) type
integer, parameter :: MNMLOB = 0, NMLOB = 0, MNMLST = 0, NMLST = 0, MNMLSX = 0, NMLSX = 0, MNMLCI = 0, NMLCI = 0, MNMLDX = 0, &
                      NMLDX = 0
integer :: NML = 0, MNML = 0, NTEST

if (type == 'OB') then
  NML = NMLOB
  MNML = MNMLOB
else if (type == 'SX') then
  NML = NMLSX
  MNML = MNMLSX
else if (type == 'DX') then
  NML = NMLDX
  MNML = MNMLDX
else if (type == 'ST') then
  NML = NMLST
  MNML = MNMLST
else if (type == 'CI') then
  NML = NMLCI
  MNML = MNMLCI
end if

if (IWAY == 1) then
  !ISM = (IPARI-1)*NML+MNML-1
  ISM = (IPARI-1)*NML+IML-MNML+1
else if (IWAY == 2) then
  if (ISM > NML) then
    IPARI = 2
    IML = ISM-NML+MNML-1
  else
    IPARI = 1
    IML = ISM+MNML-1
  end if
else
  write(6,*) ' Error in MLSM, IWAY = ',IWAY
  write(6,*) ' MLSM stop !!!'
  !stop 20
  call SYSABENDMSG('lucia_util/mlsm','Internal error','')
end if

NTEST = 0
if (NTEST /= 0) then
  write(6,'(A,A)') ' MLSM speaking, type= ',type
  write(6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
end if

end subroutine MLSM
