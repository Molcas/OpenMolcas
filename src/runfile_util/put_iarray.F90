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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine puts array integer data to the runfile.                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_iArray
!
!> @brief
!>   Add/update array data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put array data of type
!> ``Integer`` into the runfile. The data items are
!> identified by the label. Below is a list of the
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
!> - '``Center Index``'
!> - '``Ctr Index Prim``'       Idem with primitive basis set.
!> - '``nAsh``'                 The number of active orbitals per irreducible representation.
!> - '``nBas``'                 The number of basis functions per irreducible representation.
!> - '``nDel``'                 The number of deleted orbitals per irreducible representation.
!> - '``nFro``'                 The number of frozen orbitals per irreducible representation, i.e. orbitals that are not optimized.
!> - '``nIsh``'                 The number of inactive orbitals per irreducible representation.
!> - '``nIsh beta``'
!> - '``nOrb``'                 The total number of orbitals per irreducible representation.
!> - '``Orbital Type``'
!> - '``Slapaf Info 1``'        Misc. information for module SLAPAF.
!> - '``Symmetry operations``'  The symmetry operations of the point group.
!> - '``Non valence orbitals``' The total number of non valence orbitals per irreducible representation.
!> - '``MkNemo.hDisp``'         The hash matrix for displacements as specified in the mknemo module.
!>
!> @param[in] Label Name of field
!> @param[in] iData Data to put on runfile
!> @param[in] nData Length of array
!***********************************************************************

subroutine Put_iArray(Label,iData,nData)

use Definitions, only: iwp, u6

implicit none
character(len=*) :: Label
integer(kind=iwp) :: nData, iData(nData)
#include "pg_ia_info.fh"
integer(kind=iwp) :: i, item, iTmp, nTmp, RecIdx(nTocIA) = 0, RecLen(nTocIA) = 0
character(len=16) :: CmpLab1, CmpLab2, RecLab(nTocIA) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
call ffRun('iArray labels',nTmp,iTmp)
if (nTmp == 0) then
  do i=1,nTocIA
    RecLab(i) = ' '
    RecIdx(i) = sNotUsed
    RecLen(i) = 0
  end do

  ! Observe that label is at most 16 characters!

  !             1234567890123456
  RecLab(1)  = 'Center Index    '
  RecLab(2)  = 'nAsh            '
  RecLab(3)  = 'nBas            '
  RecLab(4)  = 'nDel            '
  RecLab(5)  = 'nFro            '
  RecLab(6)  = 'nIsh            '
  RecLab(7)  = 'nIsh beta       '
  RecLab(8)  = 'nOrb            '
  RecLab(9)  = 'Orbital Type    '
  RecLab(10) = 'Slapaf Info 1   '
  RecLab(11) = 'Symmetry operati' !ons
  RecLab(12) = 'nIsh_ab         '
  RecLab(13) = 'nStab           '
  RecLab(14) = '                ' !Free slot
  RecLab(15) = 'Quad_i          '
  RecLab(16) = 'RFcInfo         '
  RecLab(17) = 'RFiInfo         '
  RecLab(18) = 'RFlInfo         '
  RecLab(19) = 'SCFInfoI        '
  RecLab(20) = 'Misc            '
  RecLab(21) = 'SewIInfo        '
  RecLab(22) = '                ' !Free slot
  RecLab(23) = 'SCFInfoI_ab     '
  RecLab(24) = 'icDmp           '
  RecLab(25) = 'Symmetry Info   '
  RecLab(26) = 'Sizes           '
  RecLab(27) = '                ' !Free slot
  RecLab(28) = 'IndS            '
  RecLab(29) = '                ' !Free slot
  RecLab(30) = '                ' !Free slot
  RecLab(31) = '                ' !Free slot
  RecLab(32) = '                ' !Free Slot
  RecLab(33) = '                ' !Free slot
  RecLab(34) = '                ' !Free slot
  RecLab(35) = '                ' !Free slot
  RecLab(36) = 'LP_A            '
  RecLab(37) = 'NumCho          ' ! Number of Cholesky vectors.
  RecLab(38) = 'nFroPT          ' ! Number of Frozen for PT
  RecLab(39) = 'nDelPT          ' ! Number of Deleted for PT
  RecLab(40) = 'BasType         '
  RecLab(41) = 'Spread of Coord.'
  RecLab(42) = 'Unit Cell Atoms '
  RecLab(43) = 'iSOShl          '
  RecLab(44) = 'Non valence orbi' !tals
  RecLab(45) = 'LoProp nInts    '
  RecLab(46) = 'LoProp iSyLbl   '
  RecLab(47) = 'nDel_go         '
  RecLab(48) = 'nBas_Prim       '
  RecLab(49) = 'IsMM            '
  RecLab(50) = 'Atom -> Basis   '
  RecLab(51) = 'Logical_Info    '
  RecLab(52) = '                ' !Free slot
  RecLab(53) = '                ' !Free slot
  RecLab(54) = 'SCF nOcc        '
  RecLab(55) = 'SCF nOcc_ab     '
  RecLab(56) = '                ' !Free slot
  RecLab(57) = '                ' !Free slot
  RecLab(58) = 'iAOtSO          '
  RecLab(59) = 'iSOInf          '
  RecLab(60) = '                ' !Free Slot
  RecLab(61) = 'AuxShell        '
  RecLab(62) = 'nVec_RI         '
  RecLab(63) = 'MkNemo.hDisp    '
  RecLab(64) = 'Index ZMAT      '
  RecLab(65) = 'NAT ZMAT        '
  RecLab(66) = '                ' ! Free slot
  RecLab(67) = 'nDisp           '
  RecLab(68) = 'DegDisp         '
  RecLab(69) = 'LBList          '
  RecLab(71) = 'Ctr Index Prim  '
  RecLab(72) = 'MLTP_SINGLE     '
  RecLab(73) = 'JBNUM_SINGLE    '
  RecLab(74) = 'LROOT_SINGLE    '
  RecLab(75) = 'GeoInfo         '
  RecLab(76) = 'Cholesky BkmDim '
  RecLab(77) = 'Cholesky BkmVec '
  RecLab(78) = 'Atom Types      '
  RecLab(79) = 'LA Def          '
  RecLab(80) = 'Basis IDs       '
  RecLab(81) = 'Desym Basis IDs '
  RecLab(82) = 'primitive ids   '
  RecLab(83) = 'Root Mapping    '
  RecLab(84) = 'Fermion IDs     '
  RecLab(85) = 'IsMM Atoms      '
  RecLab(86) = 'Un_cen Charge   '
  RecLab(87) = 'PCM_N           '
  RecLab(88) = 'PCMiSph         '
  RecLab(89) = 'NVert           '
  RecLab(90) = 'IntSph          '
  RecLab(91) = 'NewSph          '
  RecLab(92) = 'XMolnr          '
  RecLab(93) = 'XEle            '
  RecLab(94) = 'iDmp            '
  RecLab(95) = 'iDmp:S          '
  Reclab(96) = 'NSTAT_SINGLE    '
  !             1234567890123456

  ! Do not go beyond 128 without changing the length of RecLab in include file too!
  call cWrRun('iArray labels',RecLab,16*nTocIA)
  call iWrRun('iArray indices',RecIdx,nTocIA)
  call iWrRun('iArray lengths',RecLen,nTocIA)
else
  call cRdRun('iArray labels',RecLab,16*nTocIA)
  call iRdRun('iArray indices',RecIdx,nTocIA)
  call iRdRun('iArray lengths',RecLen,nTocIA)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocIA
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocIA
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('iArray labels',RecLab,16*nTocIA)
    call iWrRun('iArray indices',RecIdx,nTocIA)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary iArray field'
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
if (item == -1) call SysAbendMsg('put_iArray','Could not locate',Label)
call iWrRun(RecLab(item),iData,nData)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('iArray indices',RecIdx,nTocIA)
end if
if (RecLen(item) /= nData) then
  RecLen(item) = nData
  call iWrRun('iArray lengths',RecLen,nTocIA)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Put_iArray
