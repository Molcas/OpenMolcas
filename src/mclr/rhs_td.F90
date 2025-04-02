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
use ipPage, only: ipin, W
use MCLR_Data, only: G2sq, G1t
use MCLR_Data, only: nDens, nCMO, n2Dens, ipCI, ipCM, ipMat, ipMatBA, ipMatLT, nA, nConf1, nDens2, nMBA
use MCLR_Data, only: DspVec
use MCLR_procedures, only: CISigma_td
use input_mclr, only: Debug, nSym, iMethod, State_Sym, nAsh, nBas, nIsh, nTPert
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: u6

implicit none
real*8 Temp1(nDens), Temp2(nDens), Temp3(nDens), Temp4(nDens), Temp5(nDens), Temp6(nDens), rKappa(nDens)
integer ipst, iDisp, lOper
real*8 CMO(nCMO)
integer jdisp
logical CI
character(len=8) Label
real*8 rDum(1)
real*8, allocatable :: FiX(:), MOX(:), MOT(:), MOT2(:)
integer iRC, iDSym, iOpt, iOp, ip, iS, jS, iAsh, jAsh
real*8 Dij, Ena, E2_TD

!                                                                      *
!***********************************************************************
!                                                                      *
debug = .true.
iRC = -1
idsym = loper+1
iOpt = 0
iOp = 2**loper

!----------------------------------------------------------------------*

! Read in connection matrix
! and transform it to MO basis

if (iand(ntpert(idisp),2**3) == 8) then
  iRC = -1
  iOpt = 0
  iOp = 2**loper
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
      if (ieor(iS-1,jS-1) == loper) then
        if (nBas(is)*nBas(js) /= 0) then
          if (is == js) then
            call Square(Temp6(ipMatLT(is,js)),Temp5,1,nBas(is),nBas(is))
            ip = ip+nTri_Elem(nBas(is))
          else
            call dcopy_(nBas(iS)*nBas(jS),Temp6(ipMatLt(is,js)),1,Temp5,1)
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
  call mma_allocate(MOX,n2dens,Label='MOX')
else
  call mma_allocate(MOX,1,Label='MOX')  ! Dummy allocation
end if
MOX(:) = Zero
call mma_allocate(FiX,nDens2,Label='FiX')
call IntX(FIX,Temp6,Temp5,Temp4,Temp3,rkappa,MOX,loper,idisp)

!----------------------------------------------------------------------*

!             C O N N E C T
!
!        {kappa MO}           F({D,k}){F,k}
!                                             ~
! Area for one index transformed integrals (pj|kl)

if (iand(ntpert(idisp),2**3) == 8) then
  if (iMethod == 2) then
    call mma_allocate(MOT,nmba,Label='MOT')
    call mma_allocate(MOT2,nmba,Label='MOT2')
    MOT(:) = Zero
    MOT2(:) = Zero

    ! kappa rmo Fi Fa

    ! IFG: this was outside "if (imethod == 2)",
    !      probably a bug? ipmot & ipmot2 would be uninitialized
    call r2ElInt(Temp1,MOT,MOT2,Temp3,Temp4,iDSym,One,-Half,0)
    call DaXpY_(nmba,One,MOT2,1,MOT,1)
    call mma_deallocate(MOT2)
  end if
  ! ix  ix  ~i
  call dcopy_(ndens2,[Zero],0,Temp6,1)
  ! F  =F  + F
  call DaXpY_(nDens2,One,Temp3,1,FIX,1)

  if (iMethod == 2) call CreQ_td(Temp5,MOT,G2sq,loper+1)

  do iS=1,nSym
    jS = ieor(iS-1,loper)+1
    ! F~=2*Fi~
    call DaXpY_(nIsh(is)*nBas(js),Two,Temp3(ipMat(js,is)),1,Temp6(ipMat(js,is)),1)
    if (iMethod == 2) then
      ! F~=F~+2*FA~
      call DaXpY_(nIsh(is)*nBas(js),Two,Temp4(ipMat(js,is)),1,Temp6(ipMat(js,is)),1)
      do iAsh=1,nAsh(iS)
        do jAsh=1,nAsh(is)
          Dij = G1t(iTri(iash+nA(is),jAsh+nA(is)))

          ! F~=F~+DFi~

          call DaXpY_(nBas(jS),Dij,Temp3(ipMat(js,is)+nBas(js)*(nish(is)+iAsh-1)),1, &
                      Temp6(ipMat(js,is)+nBas(js)*(nish(is)+jAsh-1)),1)
        end do
      end do
      ! F~=F~+Q~
      call DaXpY_(nAsh(is)*nBas(js),One,Temp5(ipMatba(js,is)),1,Temp6(ipMat(js,is)+nBas(js)*nIsh(is)),1)
    end if
  end do ! is

! Calculate connection contribution to hessian

end if ! ntpert
call Hess(Temp6,rkappa,Temp1,Temp3,Temp4,Temp5,Temp2,loper+1,jdisp,idisp)

! F=F~+Fx
if (iand(ntpert(idisp),2**3) == 8) call daxpy_(nDens,One,Temp6,1,rKappa,1)

! Add connection to 2el MO integrals

! Adds (pb|cd) to triangular (ab|cd)
if ((iMethod == 2) .and. (iand(ntpert(idisp),2**2) == 4)) call ABXpY(MOT,MOX,idsym)

if (CI) then
  if (iand(ntPert(idisp),2**3) /= 0) then
    call CiSigma_td(0,State_Sym,ieor(State_sym-1,idsym-1)+1,Fix,nDens2,MOX,size(MOX),rdum,1,ipCI,ipst,'N',.true.)
  else
    call CiSigma_td(0,State_Sym,ieor(State_sym-1,idsym-1)+1,Fix,nDens2,rdum,1,rdum,1,ipCI,ipst,'N',.false.)
  end if

  call ipin(ipST)
  if (idsym == 1) then
    EnA = E2_td(Fix,MOX,idsym-1,idisp)
    call ipin(ipCI)
    call DaXpY_(nConf1,-Ena,W(ipCI)%A,1,W(ipST)%A,1)
  end if
  call dscal_(nconf1,Two,W(ipST)%A,1)
end if

call DYAX(ndens2,Two,rkappa,1,Temp1,1)

do iS=1,nSym
  js = ieor(is-1,loper)+1
  if (nbas(is)*nBas(js) /= 0) &
    call DGESUB(Temp1(ipMat(is,js)),nBas(is),'N',Temp1(ipMat(js,is)),nBas(js),'T',rKappa(ipMat(is,js)),nBas(is),nBas(is),nBas(js))
end do

call mma_deallocate(FIX)
call mma_deallocate(MOX)
call mma_deallocate(MOT,safe='*')

end subroutine RHS_td
