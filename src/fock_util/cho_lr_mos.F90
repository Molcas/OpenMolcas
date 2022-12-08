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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine CHO_LR_MOs(iOK,nDen,nSym,nBas,nIsh,CM,MSQ)
!***********************************************************************
!
!   Computes Left (L) and Right (R) Cholesky vectors to represent a
!            non-symmetric density matrix
!
!            D = C[a] * C[b]^T = L * R^T
!
!   The code is generalized to treat density matrices with any number
!            (nDen) of distinct blocks. The corresponding Cholesky
!            vectors are returned in the arrays defined by the
!            DSBA_typed array MSQ
!
!   A non-zero value for the return code iOK indicates that something
!            went wrong and therefore the returned arrays may contain
!            junk
!
!   Author: F. Aquilante    (October 2008)
!
!***********************************************************************

use Data_Structures, only: DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: iOK
integer(kind=iwp), intent(in) :: nDen, nSym, nBas(nSym), nIsh(nSym)
type(DSBA_Type), intent(in) :: CM(nDen)
type(DSBA_Type), target, intent(inout) :: MSQ(nDen)
integer(kind=iwp) :: ia, ib, ikc, irc, iSym, iv, ja, jden, jv, k, l, n2b, nBm, NumV
real(kind=wp) :: Thr, Ymax
real(kind=wp), allocatable, target :: DMat0(:), PMat0(:)
real(kind=wp), pointer :: DMat(:,:) => null(), PMat(:,:,:,:) => null(), V(:,:) => null()

!***********************************************************************
irc = 0
ikc = 0

nBm = nBas(1)
do iSym=2,nSym
  nBm = max(nBm,nBas(iSym))
end do
call mma_allocate(DMat0,nBm**2,Label='DMat0')
if (nDen > 1) then
  call mma_allocate(PMat0,2*(nDen*nBm)**2,Label='PMat0')
end if

iSym = 1
do while (iSym <= nSym)

  DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)
  if (nDen == 1) then
    PMat(1:nBas(iSym),1:1,1:nBas(iSym),1:1) => DMat0(1:nBas(iSym)**2)
  else
    PMat(1:nBas(iSym),1:nDen,1:nBas(iSym),1:nDen) => PMat0(1:(nDen*nBas(iSym))**2)
  end if

  if ((nBas(iSym) > 0) .and. (nIsh(iSym) > 0)) then

    ! Inactive P[kl](a,b) = sum_i C[k](a,i)*C[l](b,i)

    n2b = nDen*nBas(iSym)

    do k=1,nDen

      do l=1,nDen

        call DGEMM_('N','T',nBas(iSym),nBas(iSym),nIsh(iSym),One,CM(k)%SB(iSym)%A2,nBas(iSym),CM(l)%SB(iSym)%A2,nBas(iSym),Zero, &
                    DMat,nBas(iSym))

        do ib=1,nBas(iSym)
          PMat(:,k,ib,l) = DMat(:,ib)
        end do

      end do

    end do

    Ymax = Zero
    do ja=1,n2b ! Loop over the compound index
      k = (ja-1)/nBas(iSym)+1
      ia = ja-(k-1)*nBas(iSym)
      Ymax = max(Ymax,PMat(ia,k,ia,k))
    end do
    Thr = 1.0e-13_wp*Ymax

    if (nDen == 1) then
      V(1:,1:) => MSQ(1)%SB(iSym)%A2(:,:)
    else
      V(1:n2b,1:n2b) => PMat0((nDen*nBm)**2+1:(nDen*nBm)**2+n2b**2)
    end if

    call CD_InCore(PMat,n2b,V,n2b,NumV,Thr,irc)

    if (NumV /= nIsh(iSym)) ikc = 1

    if (nDen > 1) then
      do jden=1,nDen
        do jv=1,NumV
          iv = 1+nBas(iSym)*(jden-1)
          call dcopy_(nBas(iSym),V(iv:,jv),1,MSQ(jDen)%SB(iSym)%A2(:,jv),1)
        end do
      end do
    end if

  end if

  if ((irc /= 0) .or. (ikc /= 0)) iSym = nSym

  iSym = iSym+1

  V => null()
  DMat => null()
  PMat => null()

end do

if (nDen > 1) call mma_deallocate(PMat0)
call mma_deallocate(DMat0)

iOK = 0
if ((irc /= 0) .or. (ikc /= 0)) iOK = 1

return

end subroutine CHO_LR_MOs
