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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine puts scalar integer data to the runfile.                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_iScalar
!
!> @brief
!>   Add/update scalar data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type
!> ``Integer`` into the runfile. The data items are
!> identified by the \p label. Below is a list of the
!> data items that are recognized. The labels are
!> case insensitive and significant to 16 characters.
!>
!> For development purposes you can use an unsupported
!> label. Whenever such a field is accessed a warning
!> message is printed in the output, to remind the
!> developer to update this routine.
!>
!> List of known labels:
!>
!> - '``Multiplicity``'               The spin multiplicity of the last SCF or RASSCF calculation.
!> - '``nMEP``'                       Number of points on the minimum energy path.
!> - '``No of Internal coordinates``' The number of internal coordinates for the molecule that is allowed within the given point
!>                                    group.
!> - '``nSym``'                       The number of irreducible representations of the molecule.
!> - '``PCM info length``'            Length of the block containing misc. info for the PCM model.
!> - '``Relax CASSCF root``'          Signals which root to perform geometry optimization for in a state average CASSCF geometry
!>                                    optimization.
!> - '``SA ready``'                   Signals that SA wavefunction is ready for gradient calculations.
!> - '``System BitSwitch``'           A bit switch controlling various functions. Will be replaced!
!> - '``Unique atoms``'
!> - '``nActel``'                     The number of active electrons in CASSCF calculation.
!> - '``MkNemo.nMole``'               The number of molecules as specified in the mknemo module.
!> - '``nLambda``'                    The number of constraints in the PCO.
!> - '``DNG``'                        Force numerical gradients.
!> - '``HessIter``'                   Last iteration where the analytical Hessian was computed.
!> - '``TS Search``'                  Flag to mark if a TS search has been activated.
!> - '``CHCCLarge``'                  Segmentation of VOs in CHCC.
!> - '``Seed``'                       The seed number for random number generator used in surface hoping.
!>
!> @param[in] Label Name of field
!> @param[in] iData Data to put on runfile
!***********************************************************************

subroutine Put_iScalar(Label,iData)

use Definitions, only: iwp, u6

implicit none
character(len=*) :: Label
integer(kind=iwp) :: iData
#include "pg_is_info.fh"
integer(kind=iwp) :: i, item, iTmp, nData, RecIdx(nTocIS) = 0, RecVal(nTocIS) = 0
character(len=16) :: CmpLab1, CmpLab2, RecLab(nTocIS) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
call ffRun('iScalar labels',nData,iTmp)
if (nData == 0) then
  do i=1,nTocIS
    RecLab(i) = ' '
    RecVal(i) = 0
    RecIdx(i) = sNotUsed
  end do

  ! Observe that label is at most 16 characters!

  !             1234567890123456
  RecLab(1)  = 'Multiplicity    '
  RecLab(2)  = 'nMEP            '
  RecLab(3)  = 'No of Internal c' !oordinates
  RecLab(4)  = 'nSym            '
  RecLab(5)  = 'PCM info length '
  RecLab(6)  = 'Relax CASSCF roo' !t
  RecLab(7)  = 'System BitSwitch'
  RecLab(8)  = 'Unique atoms    '
  RecLab(9)  = 'LP_nCenter      '
  RecLab(10) = 'ChoIni          '
  RecLab(11) = 'Unit Cell NAtoms'
  RecLab(12) = 'Cholesky Reorder'
  RecLab(13) = 'ChoVec Address  '
  RecLab(14) = 'SA ready        '
  RecLab(15) = 'NumGradRoot     '
  RecLab(16) = 'Number of roots '
  RecLab(17) = 'LoProp Restart  '
  RecLab(18) = 'MpProp nOcOb    '
  RecLab(19) = 'Highest Mltpl   '
  RecLab(20) = 'nActel          '
  RecLab(21) = 'Run_Mode        '
  RecLab(22) = 'Grad ready      '
  RecLab(23) = 'ISPIN           '
  RecLab(24) = 'SCF mode        '
  RecLab(25) = 'MkNemo.nMole    '
  RecLab(26) = 'N ZMAT          '
  RecLab(27) = 'Bfn Atoms       '
  RecLab(28) = 'FMM             '
  RecLab(29) = 'Pseudo atoms    '
  RecLab(30) = 'nChDisp         '
  RecLab(31) = 'iOff_Iter       '
  RecLab(32) = 'Columbus        '
  RecLab(33) = 'ColGradMode     '
  RecLab(34) = 'IRC             '
  RecLab(35) = 'MaxHops         '
  RecLab(36) = 'nRasHole        '
  RecLab(37) = 'nRasElec        '
  RecLab(38) = 'Rotational Symme' !try Number
  RecLab(39) = 'Saddle Iter     '
  RecLab(40) = 'iMass           '
  RecLab(41) = 'mp2prpt         ' ! True(=1) if mbpt2 was run with prpt
  RecLab(42) = 'NJOB_SINGLE     '
  RecLab(43) = 'MXJOB_SINGLE    '
  RecLab(44) = 'NSS_SINGLE      '
  RecLab(45) = 'NSTATE_SINGLE   '
  RecLab(46) = 'LDF Status      ' ! Initialized or not
  RecLab(47) = 'DF Mode         ' ! Local (1) or non-local (0) DF
  RecLab(48) = '                ' ! unused
  RecLab(49) = 'LDF Constraint  ' ! Constraint type for LDF
  RecLab(50) = 'OptimType       ' ! Optimization type in hyper
  RecLab(51) = 'STSYM           ' ! symmetry of the CAS root(s)
  RecLab(52) = 'RF CASSCF root  '
  RecLab(53) = 'RF0CASSCF root  '
  RecLab(54) = 'nCoordFiles     ' ! number of xyz-files in gateway
  RecLab(55) = 'nLambda         '
  RecLab(56) = 'DNG             '
  RecLab(57) = 'HessIter        '
  !RecLab( 58) = 'GEO_nConnect    '
  RecLab(58) = 'CHCCLarge       ' ! Segmentation of VOs in CHCC
  RecLab(59) = 'TS Search       '
  RecLab(60) = 'Number of Hops  '
  RecLab(61) = 'hopped          '
  RecLab(62) = 'Invert constrain' !ts
  RecLab(63) = 'Keep old gradien' !t
  RecLab(64) = 'embpot          ' ! Flag whether an embedding potential is present
  RecLab(65) = 'nPrim           '
  RecLab(66) = 'Seed            '
  RecLab(67) = 'Track Done      '
  RecLab(68) = 'MaxHopsTully    '
  RecLab(69) = 'EFP             ' ! Flag Effective fragment potentials
  RecLab(70) = 'nEFP_fragments  '
  RecLab(71) = 'Coor_Type       ' ! EFP fragment coordinate format
  RecLab(72) = 'nEFP_Coor       ' ! Associated number of coordinates per fragment
  RecLab(73) = 'Relax Original r' !oot
  RecLab(74) = 'Unique centers  '
  RecLab(75) = 'nXF             '
  RecLab(76) = 'CSPF            '
  RecLab(77) = 'NCONF           ' ! For MS-PDFT gradient
  RecLab(78) = 'SH RASSI run    '
  !             1234567890123456

  ! Note, when the counter here exceeds 128 update this line
  ! and the nTocIS parameter in pg_is_info.fh!

  call cWrRun('iScalar labels',RecLab,16*nTocIS)
  call iWrRun('iScalar values',RecVal,nTocIS)
  call iWrRun('iScalar indices',RecIdx,nTocIS)
else
  call cRdRun('iScalar labels',RecLab,16*nTocIS)
  call iRdRun('iScalar values',RecVal,nTocIS)
  call iRdRun('iScalar indices',RecIdx,nTocIS)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocIS
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocIS
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('iScalar labels',RecLab,16*nTocIS)
    call iWrRun('iScalar indices',RecIdx,nTocIS)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary iScalar field'
    write(u6,*) '***   Field: ',Label
    write(u6,*) '***'
#   ifndef _DEVEL_
    call AbEnd()
#   endif
  end if
end if
!----------------------------------------------------------------------*
! Write data to disk.                                                  *
!----------------------------------------------------------------------*
if (item == -1) call SysAbendMsg('put_iScalar','Could not locate',Label)
RecVal(item) = iData
call iWrRun('iScalar values',RecVal,nTocIS)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('iScalar indices',RecIdx,nTocIS)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
do i=1,num_IS_init
  if (iLbl_IS_inmem(i) == CmpLab1) then
    i_IS_inmem(i) = iData
    IS_init(i) = 1
    return
  end if
end do

return

end subroutine Put_iScalar
