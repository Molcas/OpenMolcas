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

subroutine CHO_get_MO(iOK,nDen,nSym,nBas,nIsh,CM,MSQ)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: iOK
integer(kind=iwp), intent(in) :: nDen, nSym, nBas(nSym), nIsh(nSym)
type(DSBA_Type), intent(in) :: CM(nDen)
type(DSBA_Type), intent(_OUT_) :: MSQ(nDen)
integer(kind=iwp) :: i, iComp, ikc, iOpt, irc, iSyLbl, iSym, ja, nBm, NumV
real(kind=wp) :: Thr, Ymax
character(len=8) :: Label
type(DSBA_Type) :: SMat
real(kind=wp), allocatable :: SXMat(:)
real(kind=wp), allocatable, target :: Dmat0(:)
real(kind=wp), pointer :: Dmat(:,:) => null()

!***********************************************************************
irc = 0
ikc = 0

nBm = nBas(1)
do iSym=2,nSym
  nBm = max(nBm,nBas(iSym))
end do
call mma_allocate(Dmat0,nBm**2,Label='Dmat')

iSym = 1
do while (iSym <= nSym)

  DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)

  if ((nBas(iSym) > 0) .and. (nIsh(iSym) > 0)) then

    ! Inactive D(a,b) = sum_i C(a,i)*C(b,i)

    call DGEMM_('N','T',nBas(iSym),nBas(iSym),nIsh(iSym),One,CM(1)%SB(iSYm)%A2,nBas(iSym),CM(1)%SB(iSYm)%A2,nBas(iSym),Zero,DMat, &
                nBas(iSym))

    Ymax = Zero
    do ja=1,nBas(iSym)
      Ymax = max(Ymax,DMat(ja,ja))
    end do
    Thr = 1.0e-13_wp*Ymax

    call CD_InCore(DMat,nBas(iSym),MSQ(1)%SB(iSym)%A2,nBas(iSym),NumV,Thr,irc)

    if (NumV /= nIsh(iSym)) ikc = 1

  end if

  if ((irc /= 0) .or. (ikc /= 0)) iSym = nSym

  iSym = iSym+1

  DMat => null()

end do

if ((nDen == 2) .and. (irc == 0) .and. (ikc == 0)) then

  call Allocate_DT(SMat,nBas,nBas,nSym,aCase='TRI')
  call mma_allocate(SXMat,nBm**2,Label='SXMat')

  ! Read overlap integrals (LT-storage) and get Square-storage
  iRc = -1
  iOpt = ibset(0,sNoOri)
  iComp = 1
  iSyLbl = 1
  Label = 'Mltpl  0'
  call RdOne(iRc,iOpt,Label,iComp,SMat%A0,iSyLbl)

  ! Compute  X_b[a] = C_b U_a   where  U_a = C_a^T S X_a
  ! ----------------------------------------------------
  do i=1,nSym

    DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)

    if ((nBas(i) > 0) .and. (nIsh(i) > 0)) then

      call SQUARE(SMat%SB(i)%A1,DMat,1,NBas(i),NBas(i))

      call DGEMM_('N','N',nBas(i),nIsh(i),nBas(i),One,DMat,nBas(i),MSQ(1)%SB(i)%A2,nBas(i),Zero,SXMat,nBas(i))

      DMat(1:nIsh(iSym),1:nIsh(iSym)) => DMat0(1:nIsh(iSym)**2)

      call DGEMM_('T','N',nIsh(i),nIsh(i),nBas(i),One,CM(1)%SB(i)%A2,nBas(i),SXMat,nBas(i),Zero,DMat,nIsh(i))

      !write(u6,*) ' U_a = C_a^T S X_a   for symmetry block: ',i
      !call cho_output(DMat,1,nIsh(i),1,nIsh(i),nIsh(i),nIsh(i),1,u6)

      call DGEMM_('N','N',nBas(i),nIsh(i),nIsh(i),One,CM(2)%SB(i)%A2,nBas(i),DMat,nIsh(i),Zero,MSQ(2)%SB(i)%A2,nBas(i))

    end if

    DMat => null()

  end do

  call mma_deallocate(SXMat)
  call Deallocate_DT(SMat)

end if

call mma_deallocate(DMat0)

iOK = 0
if ((irc /= 0) .or. (ikc /= 0)) iOK = 1

return

end subroutine CHO_get_MO
