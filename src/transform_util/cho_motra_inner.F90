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

subroutine Cho_MOTra_Inner(CMO,nCMOs,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BName,Do_int,ihdf5,Do_ChoInit)
! Note: frozen and deleted orbitals are not included in the transformation.

use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCMOs, nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nDel(nSym), ihdf5
real(kind=wp), intent(in) :: CMO(nCMOs)
character(len=6), intent(in) :: BName
logical(kind=iwp), intent(in) :: Do_int, Do_ChoInit
integer(kind=iwp) :: irc, iSym, iSymp, iSymq, jSym, kdisk, Lu_Xint, lXint, n, nAux(8), nOrbp, nOrbq
real(kind=wp) :: FracMem
type(DSBA_Type), target :: CMOT(1)
real(kind=wp), allocatable :: xInt(:)
integer(kind=iwp), external :: isFreeUnit
#include "chotime.fh"

n = nBas(1)**2
do iSym=2,nSym
  n = n+nBas(iSym)**2
end do
if (n /= nCMOs) then
  ! The dimension of CMO is assumed to be nBas**2 in Transp_MOs
  ! with each symmetry block starting at
  !  sym=1: 1
  !  sym=2: 1+nBas(1)**2
  !  sym=3: 1+nBas(1)**2+nBas(2)**2
  ! etc.
  ! This differs from, e.g., subroutine PriMO where each
  ! symmetry block starts at
  !  sym=1: 1
  !  sym=2: 1+nBas(1)*nOrb(1)
  !  sym=3: 1+nBas(1)*nOrb(1)+nBas(2)*nOrb(2)
  ! etc.
  ! This is only a potential problem if orbitals were removed
  ! due to linear dependence (not deleted virtual orbitals)
  call WarningMessage(2,'Cho_MOTra_: n != nCMOs')
  write(u6,*) 'n,nCMOs=',n,nCMOs
  call Abend()
end if

nAux(1:nSym) = nBas(1:nSym)-nFro(1:nSym)-nDel(1:nSym)
call Allocate_DT(CMOT(1),nAux,nBas,nSym)

call Transp_MOs(CMO,CMOT(1)%A0,nSym,nFro,nIsh,nAsh,nSsh,nBas)

timings = .true.

if (Do_int) then
  Lu_Xint = 80
  Lu_Xint = isfreeunit(Lu_Xint)
  call DaName_mf_wa(Lu_Xint,'DIAGINT')
  lXint = 0
  do jSym=1,nSym
    do iSymq=1,nSym
      nOrbq = nIsh(iSymq)+nAsh(iSymq)+nSsh(iSymq)
      iSymp = Mul(iSymq,jSym)
      if (iSymp == iSymq) then
        lXint = lXint+nOrbq*(nOrbq+1)/2
      else if (iSymp < iSymq) then
        nOrbp = nIsh(iSymp)+nAsh(iSymp)+nSsh(iSymp)
        lXint = lXint+nOrbp*nOrbq
      end if
    end do
  end do
  call mma_allocate(xInt,lXint,Label='xInt')
else
  lXint = 1
  call mma_allocate(xInt,lXint,Label='xInt')
end if

if (Do_ChoInit) then
  FracMem = Zero ! in a parallel run set it to a sensible value
  irc = 0
  call Cho_X_Init(irc,FracMem) ! initialize cholesky info
  if (irc /= 0) then
    call WarningMessage(2,'Cho_MOTra_: non-zero rc from Cho_X_Init')
    write(u6,*) 'rc=',irc
    call Abend()
  end if
end if
call CHO_TR_drv(irc,nIsh,nAsh,nSsh,CMOT,BName,Do_int,ihdf5,xInt,lXint)
if (Do_ChoInit) then
  call Cho_X_final(irc)
  if (irc /= 0) then
    call WarningMessage(2,'Cho_MOTra_: non-zero rc from Cho_X_Final')
    write(u6,*) 'rc=',irc
    call Abend()
  end if
end if

if (Do_int) then
  call GADSum(xInt,lXint)
  kdisk = 0
  call ddafile(Lu_Xint,1,Xint,lXint,kdisk)
  call daclos(Lu_Xint)
end if
call mma_deallocate(XInt)
call Deallocate_DT(CMOT(1))

return

end subroutine Cho_MOTra_Inner
