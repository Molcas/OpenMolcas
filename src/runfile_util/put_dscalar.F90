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
! This routine puts scalar double data to the runfile.                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_dScalar
!
!> @brief
!>   To add/update scalar data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type
!> ``Real*8`` into the runfile. The data items are
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
!> - '``CASDFT energy``'  Energy for the last CASDFT calculation.
!> - '``CASPT2 energy``'  Energy for the last CASPT2 calculation.
!> - '``CASSCF energy``'  Energy for the last CASSCF calculation.
!> - '``Ener_ab``'
!> - '``KSDFT energy``'   Energy for the last KS-DFT calculation.
!> - '``Last energy``'    Last energy computed.
!> - '``PC Self Energy``' Self energy for point charges.
!> - '``PotNuc``'         Nuclear repusion energy.
!> - '``RF Self Energy``' Self energy in the Kirkwood model.
!> - '``SCF energy``'     Energy for the last SCF calculation.
!> - '``EThr``'           Energy convergence threshold.
!> - '``Thrs``'
!> - '``UHF energy``'
!> - '``DFT exch coeff``' Scaling factor for exchange terms of a density
!>   functional
!> - '``DFT corr coeff``' Scaling factor for correlation terms of a
!>   density functional
!>
!> @param[in] Label Name of field
!> @param[in] rData Data to put on runfile
!***********************************************************************

subroutine Put_dScalar(Label,rData)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=*) :: Label
real(kind=wp) :: rData
#include "pg_ds_info.fh"
integer(kind=iwp) :: i, item, iTmp, nData, RecIdx(nTocDS) = 0
real(kind=wp) :: RecVal(nTocDS) = Zero
character(len=16) :: CmpLab1, CmpLab2, RecLab(nTocDS) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
! start pow mod ---
!write(6,'(3a)') 'Runfile: put_dscalar field "',Label,'"'
! end pow mod ---
call ffRun('dScalar labels',nData,iTmp)
if (nData == 0) then
  do i=1,nTocDS
    RecLab(i) = ' '
    RecVal(i) = Zero
    RecIdx(i) = sNotUsed
  end do

  ! Observe that label is at most 16 characters!

  !             1234567890123456
  RecLab(1)  = 'CASDFT energy   '
  RecLab(2)  = 'CASPT2 energy   '
  RecLab(3)  = 'CASSCF energy   '
  RecLab(4)  = 'Ener_ab         '
  RecLab(5)  = 'KSDFT energy    '
  RecLab(6)  = 'Last energy     '
  RecLab(7)  = 'PC Self Energy  '
  RecLab(8)  = 'PotNuc          '
  RecLab(9)  = 'RF Self Energy  '
  RecLab(10) = 'SCF energy      '
  RecLab(11) = 'Thrs            '
  RecLab(12) = 'UHF energy      '
  RecLab(13) = 'E_0_NN          '
  RecLab(14) = 'W_or_el         '
  RecLab(15) = 'W_or_Inf        '
  RecLab(16) = 'EThr            '
  RecLab(17) = 'Cholesky Thresho' !ld
  RecLab(18) = 'Total Nuclear Ch' !arge
  RecLab(19) = 'Numerical Gradie' !nt rDelta
  RecLab(20) = 'MpProp Energy   '
  RecLab(21) = 'UHFSPIN         '
  RecLab(22) = 'S delete thr    '
  RecLab(23) = 'T delete thr    '
  RecLab(24) = 'MD_Etot0        '
  RecLab(25) = 'MD_Time         '
  RecLab(26) = 'LDF Accuracy    '
  RecLab(27) = 'NAD dft energy  '
  RecLab(28) = 'GradLim         '
  RecLab(29) = '                ' ! Free slot
  RecLab(30) = 'Average energy  '
  RecLab(31) = 'Timestep        '
  RecLab(32) = 'MD_Etot         '
  RecLab(33) = 'Max error       '
  RecLab(34) = 'Total Charge    ' ! total number of electrons
  RecLab(35) = 'DFT exch coeff  '
  RecLab(36) = 'DFT corr coeff  '
  RecLab(37) = 'Value_l         '
  RecLab(38) = 'R_WF_HMC        '
  !             1234567890123456

  ! If you go beyond 64: update pg_ds_info.fh and this line!
  call cWrRun('dScalar labels',RecLab,16*nTocDS)
  call dWrRun('dScalar values',RecVal,nTocDS)
  call iWrRun('dScalar indices',RecIdx,nTocDS)
else
  call cRdRun('dScalar labels',RecLab,16*nTocDS)
  call dRdRun('dScalar values',RecVal,nTocDS)
  call iRdRun('dScalar indices',RecIdx,nTocDS)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocDS
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocDS
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('dScalar labels',RecLab,16*nTocDS)
    call iWrRun('dScalar indices',RecIdx,nTocDS)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary dScalar field'
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
if (item == -1) call SysAbendMsg('put_dScalar','Could not locate',Label)
RecVal(item) = rData
call dWrRun('dScalar values',RecVal,nTocDS)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('dScalar indices',RecIdx,nTocDS)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
do i=1,num_DS_init
  if (iLbl_DS_inmem(i) == CmpLab1) then
    i_DS_inmem(i) = rData
    DS_init(i) = 1
    return
  end if
end do

return

end subroutine Put_dScalar
