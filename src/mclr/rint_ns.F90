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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine RInt_ns(rkappa,rmo,Fock,Focki,idsym,reco,jspin)
!                              ~
!     Constructs  F  = <0|[E  ,H]|0>
! added in rinttd ( + <0|[[E  , Kappa],H]|0> )
!                  pq       pq                 pq
!
! Some modifications (subroutines that ends with _ns) to handle
! non anti-symmetric orbital rotations.
! Instead of using Q-Q^t as the "Q" contribution to the gradient
! Q^A-Q^{Bt} is used (see the paper). Major modifications in read2.
! Instead of constructing one set of MO integrals/particle we construct
! one set of integrals used for constructing Q^A and one for Q^B
!
! Fock is E*d/dx(lambda)
! rkappa is d/dx(lambda)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: G1t, G2sq, ipMat, ipMatBA, nA, nDens, nMBA
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rkappa(nDens), reco
real(kind=wp), intent(_OUT_) :: rMO(*)
real(kind=wp), intent(out) :: Fock(nDens), FockI(nDens)
integer(kind=iwp), intent(in) :: iDSym, jSpin
integer(kind=iwp) :: iAsh, ipF, ipFI, iS, jAsh, jS
real(kind=wp) :: Dij, Fact
real(kind=wp), allocatable :: FA(:), MT1(:), MT2(:), QA(:), QB(:)

call mma_allocate(FA,nDens,Label='FA')
! Fact controls the sign of H(k)
Fact = One
call mma_allocate(MT1,nmba,Label='MT1')
call mma_allocate(MT2,nmba,Label='MT2')
MT1(:) = Zero
MT2(:) = Zero

call R2ElInt_ns(rkappa,MT1,MT2,focki,FA,idSym,ReCo,Fact,jspin)

!if (idsym == 2) then
!  jpCMO = 1
!  do iSym=1,nSym
!    write(u6,'(A,i2.2)') 'Inactive fackmatrix = ',iSym
!    call RecPrt(' ',' ',focki(jpCMO),nBas(iSym),nBas(iSym))
!    jpCMO = jpCMO+nBas(iSym)*nBas(iSym)
!  end do
!  jpCMO = 1
!  do iSym=1,nSym
!    write(u6,'(A,i2.2)') 'Active fockmatrix     = ',iSym
!    call RecPrt(' ',' ',FA(jpCMO),nBas(iSym),nBas(iSym))
!    jpCMO = jpCMO+nBas(iSym)*nBas(iSym)
!  end do
!end if

Fock(:) = Zero

! Q  = sum(jkl)=(pj|kl)d(ijkl)
!  pi

if (iMethod == 2) then
  call mma_allocate(qA,nDens,Label='QA')
  call mma_allocate(qB,nDens,Label='QB')
  call CreQ_td(QB,MT1,G2sq,idsym)
  call CreQ_td(QA,MT2,G2sq,idsym)
end if

!call RECPRT('QB',' ',QB,nDens,1)
!call RECPRT('QA',' ',QA,nDens,1)

do iS=1,nSym
  jS = Mul(iS,idsym)

  !         I    A
  ! F  = 2 ( F  + F  )
  !  pi       pi   pi

  ! IFG This looks like a no-op. Should the second call be 'T'?
  call DGEACC(Two,Focki(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nBas(is),nIsh(js))
  call DGEACC(-Two,Focki(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nIsh(is),nBas(js))
  if (iMethod == 2) then

    ! 61-121 is all to do with multiconf

    ! IFG This looks like a no-op. Should the second call be 'T'?
    call DGEACC(Two,FA(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nBas(is),nIsh(js))
    call DGEACC(-Two,FA(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nIsh(is),nBas(js))

    do iAsh=1,nAsh(jS)
      do jAsh=1,nAsh(js)

        !         I
        ! F  = F - F  D
        !  ap   ap  ap ba

        Dij = G1t(iTri(iash+nA(js),jAsh+nA(js)))
        ipF = ipMat(is,js)+(Nish(js)+iAsh-1)*nBas(is)
        ipFI = ipMat(is,js)+(Nish(js)+jAsh-1)*nBas(is)

        !         I
        ! F  = F + F  D
        !  pa   pa  pb ab

        Fock(ipF:ipF+nBas(is)-1) = Fock(ipF:ipF+nBas(is)-1)+Dij*FockI(ipFI:ipFI+nBas(is)-1)
      end do
    end do
    do iAsh=1,nAsh(iS)
      do jAsh=1,nAsh(is)
        ipF = ipMat(is,js)+nIsh(is)+jAsh-1
        ipFI = ipMat(is,js)+nIsh(is)+iAsh-1
        Dij = G1t(iTri(iash+nA(is),jAsh+nA(is)))

        !         I
        ! F  = F - F  D
        !  pa   pa  pb ab

        Fock(ipF:ipF+nBas(js)*nBas(is)-1:nBas(is)) = Fock(ipF:ipF+nBas(js)*nBas(is)-1:nBas(is))- &
                                                     Dij*FockI(ipFI:ipFI+nBas(js)*nBas(is)-1:nBas(is))
      end do
    end do

    call DGEACC(One,QB(ipMatba(is,js)),nBas(is),'N',Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),nBas(is),nAsh(js))
    call DGESUB(Fock(ipMat(is,js)+nish(is)),nBas(is),'N',QA(ipMatba(js,is)),nBas(js),'T',Fock(ipMat(is,js)+nish(is)),nBas(is), &
                nash(is),nBas(js))
  end if
  ! Transpose ipsc2
  !call mma_allocate(T,nbas(iS)*nbas(jS),Label='T')
  !T(:) = Zero
  !call DGETMO(Fock(ipmat(is,js)),nbas(is),nbas(is),nbas(js),T,nbas(js))
  !Fock(ipmat(js,is):ipmat(js,is)+nBas(jS)*nBas(iS)-1) = T(:)
  !call mma_deallocate(T)

end do
if (imethod == 2) then
  call mma_deallocate(QA)
  call mma_deallocate(QB)
end if

Fock(:) = -Two*Fock(:)
call AddGrad(rKappa,Fock,idsym,-Two*fact)
call PickMO_td(MT1,rmo,idsym)

call mma_deallocate(MT2)
call mma_deallocate(MT1)
call mma_deallocate(FA)

end subroutine RInt_ns
