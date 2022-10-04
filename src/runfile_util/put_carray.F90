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
! This routine puts array character data to the runfile.               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_cArray
!
!> @brief
!>   Add/update array data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put array data of type
!> ``Character`` into the runfile. The data items are
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
!> '``DFT functional``'     Name of the functional used for the KS-DFT calculation.
!> '``Irreps``'             Names of the irreducible representations.
!> '``Relax Method``'       Name of the method used for geometry optimizations.
!> '``Seward Title``'       The title of the calculation as specified in module SEWARD.
!> '``Slapaf Info 3``'      Misc. information for module SLAPAF.
!> '``Unique Atom Names``'  List of the names of the symmetry unique atoms.
!> '``Unique Basis Names``' List of the basis function names.
!> '``MkNemo.lMole``'       The labels of molecules as specified in mknemo module.
!> '``MkNemo.lCluster``'    The labels of clusters as specified in mknemo module.
!> '``MkNemo.lEnergy``'     The labels of energies as specified in mknemo module.
!>
!> @param[in] Label Name of field
!> @param[in] Data  Data to put on runfile
!> @param[in] nData Length of array
!***********************************************************************

subroutine Put_cArray(Label,data,nData)

implicit none
#include "pg_ca_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
character*16 myLabel
integer nData
character*(*) data
!vv character*(*) Data(nData)
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
character*16 RecLab(nTocCA)
integer RecIdx(nTocCA)
integer RecLen(nTocCA)
save RecLab
save RecIdx
save RecLen
character*16 CmpLab1
character*16 CmpLab2
integer nTmp
integer item
integer iTmp
integer i, ilen

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
myLabel = ' '
ilen = len(Label)
ilen = min(ilen,16)
myLabel = Label(1:ilen)
call ffRun('cArray labels',nTmp,iTmp)
if (nTmp == 0) then
  do i=1,nTocCA
    RecLab(i) = ' '
    RecIdx(i) = sNotUsed
    RecLen(i) = 0
  end do

  ! Observe that label is at most 16 characters!

  !             1234567890123456
  RecLab(1)  = 'DFT functional  '
  RecLab(2)  = 'Irreps          '
  RecLab(3)  = 'Relax Method    '
  RecLab(4)  = 'Seward Title    '
  RecLab(5)  = 'Slapaf Info 3   '
  RecLab(6)  = 'Unique Atom Name' !s
  RecLab(7)  = 'Unique Basis Nam' !es
  RecLab(8)  = 'LP_L            '
  RecLab(9)  = 'MkNemo.lMole    '
  RecLab(10) = 'MkNemo.lCluster '
  RecLab(11) = 'MkNemo.lEnergy  '
  RecLab(12) = 'Symbol ZMAT     '
  RecLab(13) = 'Tinker Name     '
  RecLab(14) = 'ESPF Filename   '
  RecLab(15) = 'ChDisp          '
  RecLab(16) = 'cmass           '
  RecLab(17) = 'BirthCertificate'
  RecLab(18) = 'LastEnergyMethod'
  RecLab(19) = 'MMO Labels      '
  RecLab(20) = 'MCLR Root       '
  RecLab(21) = 'Frag_Type       ' ! EFP fragment labels
  RecLab(22) = 'ABC             ' ! EFP atom labels
  RecLab(23) = 'Un_cen Names    '
  RecLab(24) = 'cDmp            '
  RecLab(25) = 'dc: cDmp        '
  RecLab(26) = 'SymmetryCInfo   '
  RecLab(27) = 'SewardXTitle    '
  RecLab(28) = 'Align_Weights   '
  RecLab(29) = 'Quad_c          '
  !             1234567890123456
  call cWrRun('cArray labels',RecLab,16*nTocCA)
  call iWrRun('cArray indices',RecIdx,nTocCA)
  call iWrRun('cArray lengths',RecLen,nTocCA)
else
  call cRdRun('cArray labels',RecLab,16*nTocCA)
  call iRdRun('cArray indices',RecIdx,nTocCA)
  call iRdRun('cArray lengths',RecLen,nTocCA)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = myLabel
call UpCase(CmpLab1)
do i=1,nTocCA
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocCA
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = myLabel
    RecIdx(item) = sSpecialField
    call cWrRun('cArray labels',RecLab,16*nTocCA)
    call iWrRun('cArray indices',RecIdx,nTocCA)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(6,*) '***'
    write(6,*) '*** Warning, writing temporary cArray field'
    write(6,*) '***   Field: ',myLabel
    write(6,*) '***'
#   ifndef _DEVEL_
    call AbEnd()
#   endif
  end if
end if
!----------------------------------------------------------------------*
! Write data to disk.                                                  *
!----------------------------------------------------------------------*
if (item == -1) call SysAbendMsg('put_cArray','Could not locate',myLabel)
call cWrRun(RecLab(item),data,nData)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('cArray indices',RecIdx,nTocCA)
end if
if (RecLen(item) /= nData) then
  RecLen(item) = nData
  call iWrRun('cArray lengths',RecLen,nTocCA)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Put_cArray
