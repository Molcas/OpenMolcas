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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine RHS(Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,rKappa,ipst,iDisp,lOper,CMO,jdisp,CI)
!***********************************************************************
!                                                                      *
!    Purpose:                                                          *
!            Read the perturbed fock operator and one electron         *
!            hamiltonian from disk and add the connection part         *
!            to the hessian.                                           *
!                                                                      *
!     In :                                                             *
!                loper : Symmetry operator for perurbation             *
!                idisp : Perturbation component                        *
!     Out                                                              *
!                rKappa: Preconditioned RHS for the perturbation       *
!                                                                      *
!     Temporary                                                        *
!                Temp1,Temp2,Temp3...                                  *
!                                                                      *
! Author: Anders Bernhardsson, 1995                                    *
!         Theoretical Chemistry, University of Lund                    *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use ipPage, only: ipin, W
use MCLR_Data, only: DspVec, G1t, G2t, ipCI, ipCM, ipMat, ipMatBA, ipMatLT, n2Dens, nA, nCMO, nConf1, nDens, nMBA
use MCLR_procedures, only: CISigma
use input_mclr, only: Debug, iMethod, nAsh, nBas, nIsh, nOrb, nSym, nTPert, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: Temp1(nDens), Temp2(nDens), Temp3(nDens), Temp4(nDens), Temp5(nDens), Temp6(nDens), rKappa(nDens)
integer(kind=iwp), intent(in) :: ipst, iDisp, lOper, jdisp
real(kind=wp), intent(in) :: CMO(nCMO)
logical(kind=iwp), intent(in) :: CI
integer(kind=iwp) :: iAsh, iDSym, ii, ij, iOp, iOpt, ip, iRC, iS, jAsh, jS
real(kind=wp) :: Dij, E2, Ena, rDum(1)
character(len=8) :: Label
real(kind=wp), allocatable :: FIX(:), MOT(:), MOT2(:), MOX(:)

!                                                                      *
!***********************************************************************
!                                                                      *
debug = .true.
iRC = -1
idsym = loper+1
iOpt = 0
iOp = ibset(0,loper)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read in connection matrix
! and transform it to MO basis

if (btest(ntpert(idisp),3)) then
  Label = 'OVRGRD'
  call dRdMCK(iRC,iOpt,Label,DspVec(iDisp),Temp6,iop)
  if (iRc /= 0) then
    write(u6,*) 'RHS: Error reading MCKINT'
    write(u6,*) 'Label=',Label
    call Abend()
  end if

  ip = 1
  do iS=1,nSym
    do jS=1,is
      if (Mul(iS,jS) == loper+1) then
        if (nOrb(is)*nOrb(js) /= 0) then
          if (is == js) then
            call Square(Temp6(ipMatLT(is,js)),Temp5,1,nBas(is),nBas(is))
            ip = ip+nTri_Elem(nBas(is))
          else
            Temp5(1:nBas(iS)*nBas(jS)) = Temp6(ipMatLt(is,js):ipMatLt(is,js)+nBas(iS)*nBas(jS)-1)
          end if
          call DGEMM_('T','N',nOrb(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),Temp5,nBas(iS),Zero,Temp4,nOrb(iS))
          call DGEMM_('N','N',nOrb(is),nOrb(jS),nBAs(jS),One,Temp4,nOrb(iS),CMO(ipCM(jS)),nBas(jS),Zero,Temp1(ipMat(iS,jS)), &
                      nOrb(is))
          if (is /= js) then
            call DGEMM_('T','T',nOrb(jS),nBas(iS),nBAs(jS),One,CMO(ipCM(jS)),nBas(js),Temp5,nBas(iS),Zero,Temp4,nOrb(jS))
            call DGEMM_('N','N',nOrb(js),nOrb(iS),nBas(iS),One,Temp4,nOrb(jS),CMO(ipCM(iS)),nBas(iS),Zero,Temp1(ipMat(jS,iS)), &
                        nOrb(jS))
          end if

        end if
      end if
    end do
  end do
end if

!----------------------------------------------------------------------*

! Read in derivative of hamiltonian

if ((iMethod == 2) .and. (n2Dens /= 0)) then
  call mma_allocate(MOX,n2Dens,Label='MOX')
else
  call mma_allocate(MOX,1,Label='MOX')
end if
MOX(:) = Zero
call mma_allocate(FiX,nDens,Label='FIX')

call IntX(FIX,Temp6,Temp5,Temp4,Temp3,rkappa,MOX,loper,idisp)

!----------------------------------------------------------------------*

!             C O N N E C T
!
!        {kappa MO}           F({D,k}){F,k}
!                                             ~
! Area for one index transformed integrals (pj|kl)

if (btest(ntpert(idisp),3)) then
  if (iMethod == 2) then
    call mma_allocate(MOT,nmba,Label='MOT')
    call mma_allocate(MOT2,nmba,Label='MOT2')
  else
    call mma_allocate(MOT,1,Label='MOT')
    call mma_allocate(MOT2,1,Label='MOT2')
  end if
  MOT(:) = Zero
  MOT2(:) = Zero

  ! kappa rmo Fi Fa

  call r2ElInt(Temp1,MOT,MOT2,Temp3,Temp4,iDSym,One,-Half,0)

  if (imethod == 2) MOT(:) = MOT(:)+MOT2(:)
  call mma_deallocate(MOT2)

  ! ix  ix  ~i
  Temp6(:) = Zero
  ! F  =F  + F
  FIX(:) = FIX(:)+Temp3(:)

  if (iMethod == 2) call CreQ(Temp5,MOT,G2t,loper+1)

  do iS=1,nSym
    jS = Mul(iS,loper+1)
    if (nOrb(js) < 1) cycle
    ! F~=2*Fi~
    Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1) = Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1)+ &
                                                           Two*Temp3(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1)
    if (iMethod == 2) then
      ! F~=F~+2*FA~
      Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1) = Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1)+ &
                                                             Two*Temp4(ipMat(js,is):ipMat(js,is)+nIsh(is)*nOrb(js)-1)
      if (nAsh(iS) < 1) cycle
      do iAsh=1,nAsh(iS)
        ii = ipMat(js,is)+nOrb(js)*(nish(is)+iAsh-1)
        do jAsh=1,nAsh(is)
          ij = ipMat(js,is)+nOrb(js)*(nish(is)+jAsh-1)
          Dij = G1t(iTri(iash+nA(is),jAsh+nA(is)))

          ! F~=F~+DFi~

          Temp6(ij:ij+nOrb(jS)-1) = Temp6(ij:ij+nOrb(jS)-1)+Dij*Temp3(ii:ii+nOrb(jS)-1)
        end do
      end do
      ! F~=F~+Q~
      Temp6(ipMat(js,is)+nOrb(js)*nIsh(is):ipMat(js,is)+nOrb(js)*(nIsh(is)+nAsh(is))-1) = &
        Temp6(ipMat(js,is)+nOrb(js)*nIsh(is):ipMat(js,is)+nOrb(js)*(nIsh(is)+nAsh(is))-1)+ &
        Temp5(ipMatba(js,is):ipMatba(js,is)+nAsh(is)*nOrb(js)-1)
    end if
  end do ! is

! Calculate connection contribution to hessian

end if ! ntpert

call Hess(Temp6,rkappa,Temp1,Temp3,Temp4,Temp5,Temp2,loper+1,jdisp,idisp)

! F=F~+Fx
if (btest(ntpert(idisp),3)) rKappa(:) = rKappa(:)+Temp6(:)

! Add connection to 2el MO integrals

! Adds (pb|cd) to triangular (ab|cd)
if ((iMethod == 2) .and. btest(ntpert(idisp),2)) call ABXpY(MOT,MOX,idsym)

if (CI) then

  call CiSigma(0,State_Sym,Mul(State_sym,idsym),FIX,nDens,MOX,size(MOX),rdum,1,ipCI,ipst,.true.)

  call ipin(ipst)
  if (idsym == 1) then
    EnA = E2(Fix,MOX,idsym-1,idisp)
    call ipin(ipCI)
    W(ipST)%A(1:nConf1) = W(ipST)%A(1:nConf1)-Ena*W(ipCI)%A(1:nConf1)
  end if
  W(ipST)%A(1:nConf1) = Two*W(ipST)%A(1:nConf1)
end if

Temp1(:) = Two*rKappa(:)

do iS=1,nSym
  js = Mul(is,loper+1)
  if (nOrb(is)*nOrb(js) /= 0) &
    call DGESUB(Temp1(ipMat(is,js)),nOrb(is),'N',Temp1(ipMat(js,is)),nOrb(js),'T',rKappa(ipMat(is,js)),nOrb(is),nOrb(is),nOrb(js))
end do

call mma_deallocate(FIX)
call mma_deallocate(MOX,safe='*')
call mma_deallocate(MOT,safe='*')
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine RHS
