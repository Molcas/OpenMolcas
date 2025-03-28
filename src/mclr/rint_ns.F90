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

use Arrays, only: G2sq, G1t
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use MCLR_Data, only: nDens2, ipMat, ipMatBA, nA, nMBA
use input_mclr, only: iMethod, nSym, nAsh, nBas, nIsh

implicit none
real*8 rkappa(nDens2), rMO(*), Fock(nDens2), FockI(ndens2)
integer iDSym, jSpin
real*8 reco
real*8, allocatable :: FA(:), MT1(:), MT2(:), QA(:), QB(:)
real*8 Fact, Dij
integer iS, jS, iAsh, jAsh, ipF, ipFI
! Statement function
integer i, j, itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

call mma_allocate(FA,ndens2,Label='FA')
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
  call mma_allocate(qA,ndens2,Label='QA')
  call mma_allocate(qB,ndens2,Label='QB')
  call CreQ_td(QB,MT1,G2sq,idsym)
  call CreQ_td(QA,MT2,G2sq,idsym)
end if

!call RECPRT('QB',' ',QB,nDens2,1)
!call RECPRT('QA',' ',QA,nDens2,1)

do iS=1,nSym
  jS = ieor(iS-1,idsym-1)+1

  !         I    A
  ! F  = 2 ( F  + F  )
  !  pi       pi   pi

  call DGEADD2(Two,Focki(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nbas(is),nish(js))
  call DGEADD2(-Two,Focki(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nish(is),nBas(js))
  if (iMethod == 2) then

    ! 61-121 is all to do with multiconf

    call DGEADD2(Two,FA(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nbas(is),nIsh(js))
    call DGEADD2(-Two,FA(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),'N',Fock(ipMat(is,js)),nBas(is),nish(is),nBas(js))

    do iAsh=1,nAsh(jS)
      do jAsh=1,nAsh(js)

        !         I
        ! F  = F - F  D
        !  ap   ap  ap ba

        Dij = G1t(itri(iash+nA(js),jAsh+nA(js)))
        ipF = ipMat(is,js)+(Nish(js)+iAsh-1)*nBas(is)
        ipFI = ipMat(is,js)+(Nish(js)+jAsh-1)*nBas(is)

        !         I
        ! F  = F + F  D
        !  pa   pa  pb ab

        call DaXpY_(nBas(is),Dij,focki(ipFI),1,Fock(ipF),1)
      end do
    end do
    do iAsh=1,nAsh(iS)
      do jAsh=1,nAsh(is)
        ipF = ipMat(is,js)+nIsh(is)+jAsh-1
        ipFI = ipMat(is,js)+nIsh(is)+iAsh-1
        Dij = G1t(itri(iash+nA(is),jAsh+nA(is)))

        !         I
        ! F  = F - F  D
        !  pa   pa  pb ab

        call DaXpY_(nBas(js),-Dij,focki(ipFI),nbas(is),Fock(ipF),nbas(is))
      end do
    end do

    call DGEACC(Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),QB(ipMatba(is,js)),nBas(is),nBas(is),nAsh(js))
    call DGESUB(Fock(ipMat(is,js)+nish(is)),nBas(is),'N',QA(ipMatba(js,is)),nBas(js),'T',Fock(ipMat(is,js)+nish(is)),nBas(is), &
                nash(is),nBas(js))
  end if
  ! Transpose ipsc2
  !call mma_allocate(T,nbas(is)*nbas(jS),Label='T')
  !call dcopy_(nbas(is)*nbas(jS),[Zero],0,T,1)
  !call DGETMO(Fock(ipmat(is,js)),nbas(is),nbas(is),nbas(js),T,nbas(js))
  !call dcopy_(nBas(jS)*nBas(iS),T,1,Fock(ipmat(js,is)),1)
  !call mma_deallocate(T)

end do
if (imethod == 2) then
  call mma_deallocate(QA)
  call mma_deallocate(QB)
end if

call DSCAL_(ndens2,-Two,fock,1)
call AddGrad(rKappa,Fock,idsym,Two*(-fact))
call PickMO_td(MT1,rmo,idsym)

call mma_deallocate(MT2)
call mma_deallocate(MT1)
call mma_deallocate(FA)

end subroutine RInt_ns
