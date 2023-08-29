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

subroutine CHO_MCA_INIT(SKIP_PRESCREEN)
!
! Purpose: initialization of Cholesky decomposition in MOLCAS.

use Index_Functions, only: nTri_Elem
use index_arrays, only: iSO2Sh
use Cholesky, only: iBas, iBasSh, IFCSew, iShlSO, iShP2Q, iShP2RS, iSOShl, iSP2F, LuPri, Mx2Sh, MxOrSh, nBas, nBasSh, nBasT, &
                    nBstSh, nnShl, nnShl_Tot, nShell, nSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: SKIP_PRESCREEN
integer(kind=iwp) :: I, IA, IJ, IJSHL, ISHL, ISYM, J, NUMIJ
character(len=*), parameter :: SECNAM = 'CHO_MCA_INIT'

! Check that the number of shells is within limits.
! -------------------------------------------------

if (NSHELL < 1) then
  write(LUPRI,*) 'NSHELL out of bounds: ',NSHELL
  call CHO_QUIT('NSHELL out of bounds in '//SECNAM,102)
end if

! Compute total #shell pair.
! --------------------------

NNSHL_TOT = nTri_Elem(NSHELL)
if (NNSHL_TOT < 1) then
  write(LUPRI,*) 'NNSHL_TOT=nTri_Elem(NSHELL) is non-positive: ',NNSHL_TOT
  write(LUPRI,*) 'Integer overflow ?'
  call CHO_QUIT('NNSHL_TOT out of bounds in '//SECNAM,102)
end if

! Compute contributing #shell pair by prescreening (if requested).
! iSP2F(k): returns global shell pair of contributing shell pair k.
!           (Allocated and defined in CHO_DIAPS.)
! -----------------------------------------------------------------

if (SKIP_PRESCREEN) then
  if ((NNSHL < 1) .or. (NNSHL > NNSHL_TOT)) then
    write(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',SKIP_PRESCREEN
    write(LUPRI,*) 'NNSHL is out-of-bounds: ',NNSHL
    write(LUPRI,*) 'Condition: 0 < NNSHL < ',NNSHL_TOT
    call CHO_QUIT('Initialization error in '//SECNAM,102)
  end if
  if (size(iSP2F) /= NNSHL) then
    write(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',SKIP_PRESCREEN
    write(LUPRI,*) 'NNSHL is: ',NNSHL
    write(LUPRI,*) 'SIZE(iSP2F) must be equal to NNSHL, SIZE(iSP2F)= ',size(iSP2F)
    call CHO_QUIT('Initialization error in '//SECNAM,102)
  end if
else
  call CHO_DIASP()
end if

! Get the number of symmetries.
! -----------------------------

call GET_ISCALAR('nSym',NSYM)  ! Get # irreps.
if ((NSYM < 1) .or. (NSYM > 8)) then
  write(LUPRI,*) 'NSYM out of bounds: ',NSYM
  call CHO_QUIT('NSYM out of bounds in '//SECNAM,102)
end if

! NBAS(ISYM): # basis functions (SOs) in symmetry ISYM
! IBAS(ISYM): offset to basis functions in symmetry ISYM
! NBAST     : total number of basis functions
! ------------------------------------------------------

call GET_IARRAY('nBas',NBAS,NSYM)
IBAS(1) = 0
NBAST = NBAS(1)
do ISYM=2,NSYM
  IBAS(ISYM) = NBAST
  NBAST = NBAST+NBAS(ISYM)
end do
if (NBAST < 1) then
  write(LUPRI,*) 'NBAST out of bounds: ',NBAST
  call CHO_QUIT('NBAST out of bounds in '//SECNAM,102)
end if

! Allocate shell based index arrays.
! ----------------------------------

call mma_allocate(iBasSh,nSym,nShell,Label='iBasSh')
call mma_allocate(nBasSh,nSym,nShell,Label='nBasSh')
call mma_allocate(nBstSh,nShell,Label='nBstSh')

! ISOSHL(I): shell to which SO I belongs
! --------------------------------------

call mma_allocate(iSOShl,NBAST,Label='iSOShl')
do ISYM=1,NSYM
  do IA=1,NBAS(ISYM)
    I = IBAS(ISYM)+IA
    iSOShl(I) = ISO2SH(I)
  end do
end do

! NBASSH(ISYM,ISHL): total dimension of shell ISHL, sym. ISYM
! NBSTSH(ISHL): total dimension of shell ISHL
! MXORSH      : max. shell dimension
! -----------------------------------------------------------

call CHO_SETSH(IBASSH,NBASSH,NBSTSH,IBAS,NBAS,ISOSHL,NSYM,NSHELL,NBAST)

MXORSH = NBSTSH(1)
do ISHL=2,NSHELL
  MXORSH = max(MXORSH,NBSTSH(ISHL))
end do

! MX2SH: max. dimension of contributing shell pair.
! -------------------------------------------------

MX2SH = -1
do IJSHL=1,NNSHL
  IJ = iSP2F(IJSHL)
  call CHO_INVPCK(IJ,I,J,.true.)
  if (I == J) then
    NUMIJ = nTri_Elem(NBSTSH(I))
  else
    NUMIJ = NBSTSH(I)*NBSTSH(J)
  end if
  MX2SH = max(MX2SH,NUMIJ)
end do
if (MX2SH < 1) then
  write(LUPRI,*) 'Max. shell pair dimension non-positive: ',MX2SH
  call CHO_QUIT('Initialization problem in '//SECNAM,102)
end if

! If needed, allocate memory for extracting qualified columns
! directly in reduced set from Seward.
! -----------------------------------------------------------

if (IFCSEW == 2) then
  call mma_allocate(iShP2RS,2,Mx2Sh,Label='iShP2RS')
  call mma_allocate(iShP2Q,2,Mx2Sh,Label='iShP2Q ')
end if

! ISHLSO(I): index of SO I within its shell
! -----------------------------------------

call mma_allocate(iShlSO,nBasT,Label='iShlSO')
call CHO_SETSH2(iShlSO,iSOShl,NBSTSH,NBAST,NSHELL)

end subroutine CHO_MCA_INIT
