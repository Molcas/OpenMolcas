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

subroutine RHS_td(Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,rKappa,ipst,iDisp,lOper,CMO,jdisp,CI)
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
use MCLR_Data, only: DspVec, G1t, G2sq, ipCI, ipCM, ipMat, ipMatBA, ipMatLT, n2Dens, nA, nCMO, nConf1, nDens, nMBA
use MCLR_procedures, only: CISigma_td
use input_mclr, only: Debug, iMethod, nAsh, nBas, nIsh, nSym, nTPert, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: Temp1(nDens), Temp2(nDens), Temp3(nDens), Temp4(nDens), Temp5(nDens), Temp6(nDens), rKappa(nDens)
integer(kind=iwp), intent(in) :: ipst, iDisp, lOper, jdisp
real(kind=wp), intent(in) :: CMO(nCMO)
logical(kind=iwp), intent(in) :: CI
integer(kind=iwp) :: iAsh, iDSym, ii, ij, iOp, iOpt, ip, iRC, iS, jAsh, jS
real(kind=wp) :: Dij, E2_TD, Ena, rDum(1)
character(len=8) :: Label
real(kind=wp), allocatable :: FiX(:), MOT(:), MOT2(:), MOX(:)

!                                                                      *
!***********************************************************************
!                                                                      *
debug = .true.
iRC = -1
idsym = loper+1
iOpt = 0
iOp = ibset(0,loper)

!----------------------------------------------------------------------*

! Read in connection matrix
! and transform it to MO basis

if (btest(ntpert(idisp),3)) then
  Label = 'OVRGRD'
  call dRdMCK(iRC,iOpt,Label,DspVec(iDisp),Temp6,iop)
  if (iRc /= 0) then
    write(u6,*)
    write(u6,*) ' *** Error in subroutine RHS_TD ***'
    write(u6,*) ' Error when reading OVRGRD from MCKINT'
    write(u6,*)
    return
  end if
  ip = 1
  do iS=1,nSym
    do jS=1,is
      if (Mul(iS,jS) == loper+1) then
        if (nBas(is)*nBas(js) /= 0) then
          if (is == js) then
            call Square(Temp6(ipMatLT(is,js)),Temp5,1,nBas(is),nBas(is))
            ip = ip+nTri_Elem(nBas(is))
          else
            Temp5(1:nBas(iS)*nBas(jS)) = Temp6(ipMatLt(is,js):ipMatLt(is,js)+nBas(iS)*nBas(jS)-1)
          end if
          call DGEMM_('T','N',nBas(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),Temp5,nBas(iS),Zero,Temp4,nBas(iS))
          call DGEMM_('N','N',nBas(is),nBas(jS),nBAs(jS),One,Temp4,nBas(iS),CMO(ipCM(jS)),nBas(jS),Zero,Temp1(ipMat(iS,jS)), &
                      nBas(is))
          if (is /= js) then
            call DGEMM_('T','T',nBas(jS),nBas(iS),nBAs(jS),One,CMO(ipCM(jS)),nBas(js),Temp5,nBas(iS),Zero,Temp4,nBas(jS))
            call DGEMM_('N','N',nBas(js),nBas(iS),nBas(iS),One,Temp4,nBas(jS),CMO(ipCM(iS)),nBas(iS),Zero,Temp1(ipMat(jS,iS)), &
                        nBas(jS))
          end if

        end if
      end if
    end do
  end do
end if

!----------------------------------------------------------------------*

! Read in derivative of hamiltonian

if (iMethod == 2) then
  call mma_allocate(MOX,n2Dens,Label='MOX')
else
  call mma_allocate(MOX,1,Label='MOX')  ! Dummy allocation
end if
MOX(:) = Zero
call mma_allocate(FiX,nDens,Label='FiX')
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
    MOT(:) = Zero
    MOT2(:) = Zero

    ! kappa rmo Fi Fa

    ! IFG: this was outside "if (imethod == 2)",
    !      probably a bug? ipmot & ipmot2 would be uninitialized
    call r2ElInt(Temp1,MOT,MOT2,Temp3,Temp4,iDSym,One,-Half,0)
    MOT(:) = MOT(:)+MOT2(:)
    call mma_deallocate(MOT2)
  end if
  ! ix  ix  ~i
  Temp6(:) = Zero
  ! F  =F  + F
  FIX(:) = FIX(:)+Temp3(:)

  if (iMethod == 2) call CreQ_td(Temp5,MOT,G2sq,loper+1)

  do iS=1,nSym
    jS = Mul(iS,loper+1)
    ! F~=2*Fi~
    Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1) = Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1)+ &
                                                           Two*Temp3(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1)
    if (iMethod == 2) then
      ! F~=F~+2*FA~
      Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1) = Temp6(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1)+ &
                                                             Two*Temp4(ipMat(js,is):ipMat(js,is)+nIsh(is)*nBas(js)-1)
      do iAsh=1,nAsh(iS)
        do jAsh=1,nAsh(is)
          Dij = G1t(iTri(iash+nA(is),jAsh+nA(is)))

          ! F~=F~+DFi~

          ii = ipMat(js,is)+nBas(js)*(nish(is)+iAsh-1)
          ij = ipMat(js,is)+nBas(js)*(nish(is)+jAsh-1)
          Temp6(ij:ij+nBas(jS)-1) = Temp6(ij:ij+nBas(jS)-1)+Dij*Temp3(ii:ii+nBas(jS)-1)
        end do
      end do
      ! F~=F~+Q~
      Temp6(ipMat(js,is)+nBas(js)*nIsh(is):ipMat(js,is)+nBas(js)*(nIsh(is)+nAsh(is))-1) = &
        Temp6(ipMat(js,is)+nBas(js)*nIsh(is):ipMat(js,is)+nBas(js)*(nIsh(is)+nAsh(is))-1)+ &
        Temp5(ipMatba(js,is):ipMatba(js,is)+nAsh(is)*nBas(js)-1)
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
  if (btest(ntPert(idisp),3)) then
    call CiSigma_td(0,State_Sym,Mul(State_sym,idsym),Fix,nDens,MOX,size(MOX),rdum,1,ipCI,ipst,'N',.true.)
  else
    call CiSigma_td(0,State_Sym,Mul(State_sym,idsym),Fix,nDens,rdum,1,rdum,1,ipCI,ipst,'N',.false.)
  end if

  call ipin(ipST)
  if (idsym == 1) then
    EnA = E2_td(Fix,MOX,idsym-1,idisp)
    call ipin(ipCI)
    W(ipST)%A(1:nConf1) = W(ipST)%A(1:nConf1)-Ena*W(ipCI)%A(1:nConf1)
  end if
  W(ipST)%A(1:nConf1) = Two*W(ipST)%A(1:nConf1)
end if

Temp1(:) = Two*rKappa(:)

do iS=1,nSym
  js = Mul(is,loper+1)
  if (nbas(is)*nBas(js) /= 0) &
    call DGESUB(Temp1(ipMat(is,js)),nBas(is),'N',Temp1(ipMat(js,is)),nBas(js),'T',rKappa(ipMat(is,js)),nBas(is),nBas(is),nBas(js))
end do

call mma_deallocate(FIX)
call mma_deallocate(MOX)
call mma_deallocate(MOT,safe='*')

end subroutine RHS_td
