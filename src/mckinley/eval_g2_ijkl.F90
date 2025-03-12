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
! Copyright (C) 1990, Roland Lindh                                     *
!               1995,1996, Anders Bernhardsson                         *
!***********************************************************************

subroutine Eval_g2_ijkl(iS,jS,kS,lS,Hess,nHess,Post_Process,iInt,n_Int,nACO,lGrad,lHess,lPick,nBuffer,Buffer,nDens,DTemp,DInAc, &
                        MOip,nTwo2,nQuad)

use setup, only: nAux
use McKinley_global, only: nMethod, RASSCF
use Index_Functions, only: iTri
use iSD_data, only: iSD, nSD
use k2_arrays, only: Aux, Create_Braket, Destroy_Braket, nFT, Sew_Scr
use stdalloc, only: mma_allocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, nHess, n_Int, nACO, nBuffer, nDens, MOip(0:7), nTwo2, nQuad
real(kind=wp), intent(inout) :: Hess(nHess), iInt(n_Int), Buffer(nBuffer), DTemp(nDens), DInAc(nDens)
logical(kind=iwp), intent(inout) :: Post_Process
logical(kind=iwp), intent(in) :: lGrad, lHess, lPick
integer(kind=iwp) :: iAng(4), iBasAO, iBasi, iBasn, iBsInc, iDer, iFnc(4), ipFin, ipMem2, ipMem3, ipMem4, ipMemX, ipMOC, ipPSO, &
                     iSD4(0:nSD,4), jBasAO, jBasj, jBasn, jBsInc, JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), kBasAO, kBask, kBasn, &
                     kBsInc, kCmp, lBasAO, lBasl, lBasn, lBsInc, lCmp, Mem1, Mem2, Mem3, Mem4, MemCMO, MemFck, MemFin, MemMax, &
                     MemPrm, MemPSO, MemX, nijkl, nRys, nSO, nTemp
real(kind=wp) :: Coor(3,4), PMax
logical(kind=iwp) :: JfG(4), JfGrd(3,4), JfHss(4,3,4,3), lDot, lDot2, lTri
real(kind=wp), pointer :: Fin(:), MOC(:), PSO(:,:), Temp(:), Work2(:), Work3(:), Work4(:), WorkX(:)
logical(kind=iwp), parameter :: n8 = .true.
integer(kind=iwp), external :: MemSO2_P

iFnc(:) = -99
PMax = Zero
if (.not. allocated(Sew_Scr)) then
  call mma_MaxDBLE(MemMax)
  if (MemMax > 8000) MemMax = MemMax-8000
  call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
else
  MemMax = size(Sew_Scr)
end if
ipMOC = 1

call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

call Coor_setup(iSD4,nSD,Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! The code is working in such away that the MO needs upper and lower
! triangular parts of ij kl but hessian needs only lower, check if the
! integralbatch is lower or upper!!

lTri = (iTri(iS,jS) >= iTri(kS,lS))
if ((.not. lTri) .and. (nMethod /= RASSCF)) return
lDot = (lTri .and. lHess)

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for zeta, eta, kappa, P and Q.
! Allocate also for Alpha, Beta , Gamma and Delta in expanded form.

call Create_BraKet(iSD4(5,1)*iSD4(5,2),iSD4(5,3)*iSD4(5,4))
!                                                                      *
!***********************************************************************
!                                                                      *
! Fix the 1st order density matrix

! Pick up pointers to desymmetrized 1st order density matrices.
! Observe that the desymmetrized 1st order density matrices
! follow the contraction index.

if (lTri .and. lPick) call Dens_Infos(nMethod)

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute total size of the second order density matrix in SO basis.
!----------------------------------------------------------------------*
nSO = MemSO2_P(nSD,iSD4)
ldot2 = ldot
if (nSO == 0) ldot2 = .false.

! Compute memory request for the primitives.

iDer = 2
if (.not. ldot2) iDer = 1
iAng(:) = iSD4(1,:)
call MemRg2(iAng,nRys,MemPrm,iDer)

!----------------------------------------------------------------------*
!
! Calculate which derivatives should be made.
!
!----------------------------------------------------------------------*

call DerCtr(ldot2,JfGrd,JndGrd,JfHss,JndHss,JfG,nSD,iSD4)

!----------------------------------------------------------------------*
!
! Decide on the partitioning of the shells based on the
! available memory and the requested memory.
!
!----------------------------------------------------------------------*

call PSOAO2(nSO,MemPrm,MemMax,iFnc,nAco,Mem1,Mem2,Mem3,Mem4,MemX,MemPSO,MemFck,nFT,MemFin,nBuffer,nSD,iSD4)

iBasi = iSD4(3,1)
jBasj = iSD4(3,2)
kBask = iSD4(3,3)
lBasl = iSD4(3,4)

iBsInc = iSD4(4,1)
jBsInc = iSD4(4,2)
kBsInc = iSD4(4,3)
lBsInc = iSD4(4,4)

kCmp = iSD(2,kS)
lCmp = iSD(2,lS)

!----------------------------------------------------------------------*
!
! Loop over basis function if we do not have enough of memory to
! calculate them in one step.
!
!----------------------------------------------------------------------*
do iBasAO=1,iBasi,iBsInc
  iBasn = min(iBsInc,iBasi-iBasAO+1)
  iSD4(8,1) = iBasAO-1
  iSD4(19,1) = iBasn

  !--------------------------------------------------------------------*
  !
  ! Move appropriate portions of the desymmetrized 1st order density matrix.
  !
  !--------------------------------------------------------------------*
  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iSD4(8,2) = jBasAO-1
    iSD4(19,2) = jBasn

    if (lpick) call Picky_Mck(nSD,iSD4,1,2,nMethod)

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iSD4(8,3) = kBasAO-1
      iSD4(19,3) = kBasn

      if (lpick) then
        call Picky_Mck(nSD,iSD4,1,3,nMethod)
        call Picky_Mck(nSD,iSD4,2,3,nMethod)
      end if

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iSD4(8,4) = lBasAO-1
        iSD4(19,4) = lBasn

        if (lpick) then
          call Picky_Mck(nSD,iSD4,3,4,nMethod)
          call Picky_Mck(nSD,iSD4,1,4,nMethod)
          call Picky_Mck(nSD,iSD4,2,4,nMethod)
        end if

        !--------------------------------------------------------------*

        nijkl = iBasn*jBasn*kBasn*lBasn

        ! Mark out the memory allocations explicitly with pointers
        ! MO tranformation buffer
        MemCMO = nACO*(kCmp*kBasn+lCmp*lBasn)
        MOC(1:MemCMO) => Sew_Scr(ipMOC:ipMOC+MemCMO-1)
        ! Area for the AO integrals
        ipFin = ipMOC+MemCMO
        Fin(1:MemFin) => Sew_Scr(ipFin:ipFin+MemFin-1)
        ! Area for 2el density
        ipPSO = ipFin+MemFin
        PSO(1:nijkl,1:nSO) => Sew_Scr(ipPSO:ipPSO+nijkl*nSO-1)
        if (nijkl*nSO > Mem1) then
          write(u6,'(A)') 'nijkl*nSO>Mem1'
          write(u6,*) 'njikl,nSO=',nijkl,nSO
          write(u6,*) 'Mem1=',Mem1
          call Abend()
        end if
        ipMem2 = ipPSO+Mem1  ! Work
        Work2(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem2-1)
        ipMem3 = ipMem2+Mem2 ! Work
        Work3(1:Mem3) => Sew_Scr(ipMem3:ipMem3+Mem3-1)
        ipMemX = ipMem3+Mem3 ! Work
        WorkX(1:MemX) => Sew_Scr(ipMemX:ipMemX+MemX-1)

        ! If MO transformation is performed in the standard way
        ! reserve memory for partial transformed integrals

        ! Multilayer

        ipMem4 = ipMem2+Mem2-Mem4
        Work4(1:Mem4) => Sew_Scr(ipMem4:ipMem4+Mem4-1)
        nTemp = Mem2+Mem3+MemX
        Temp(1:nTemp) => Sew_Scr(ipMem2:ipMem2+nTemp-1)

        !--------------------------------------------------------------*
        !
        ! Get the 2nd order density matrix in SO basis.
        !
        !--------------------------------------------------------------*

        if (n8) call PickMO(MOC,MemCMO,nSD,iSD4)
        if (ldot2) call PGet0(nijkl,PSO,nSO,iFnc,MemPSO,Temp,nTemp,nQuad,PMax,iSD4)

        ! Compute gradients of shell quadruplet

        call TwoEl_mck(Coor,nRys,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,JfG,PSO,nijkl,nSO,Work2,Mem2,Work3,Mem3,Work4,Mem4,Aux,nAux, &
                       WorkX,MemX,Fin,MemFin,Temp,nTemp,nTwo2,nFT,iInt,Buffer,nBuffer,lgrad,ldot2,n8,ltri,DTemp,DInAc,moip,nAco, &
                       MOC,MemCMO,iSD4)
        Post_Process = .true.

        nullify(MOC)
        nullify(Fin)
        nullify(PSO)
        nullify(Work2)
        nullify(Work3)
        nullify(WorkX)
        nullify(Work4)
        nullify(Temp)

        !--------------------------------------------------------------*

      end do
    end do
  end do
end do
call Destroy_Braket()

contains

subroutine Dens_Infos(nMethod)

  use Dens_stuff, only: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, &
                        ipDij, ipDij2, ipDik, ipDik2, ipDil, ipDil2, ipDjk, ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, mDCRij, mDCRik, &
                        mDCRil, mDCRjk, mDCRjl, mDCRkl
  use k2_arrays, only: ipDijS, ipDijS2
  use Index_Functions, only: iTri

  integer(kind=iwp), intent(in) :: nMethod
  integer(kind=iwp) :: ijS, ikS, ilS, ipDum, ipTmp, ipTmp2, iS, jkS, jlS, jS, klS, kS, lS
  integer(kind=iwp), parameter :: Nr_of_Densities = 1

  iS = iSD4(11,1)
  jS = iSD4(11,2)
  kS = iSD4(11,3)
  lS = iSD4(11,4)

  ijS = iTri(iS,jS)
  klS = iTri(kS,lS)
  ikS = iTri(iS,kS)
  ilS = iTri(iS,lS)
  jkS = iTri(jS,kS)
  jlS = iTri(jS,lS)
  ijS = iTri(iS,jS)
  klS = iTri(kS,lS)
  ikS = iTri(iS,kS)
  ilS = iTri(iS,lS)
  jkS = iTri(jS,kS)
  jlS = iTri(jS,lS)
  ipTmp = ipDijs
  if (nMethod == RASSCF) ipTmp2 = ipDijs2
  call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDij2,ipDDij2)
  call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDkl2,ipDDkl2)
  call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDik2,ipDDik2)
  call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDil2,ipDDil2)
  call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDjk2,ipDDjk2)
  call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDjl2,ipDDjl2)

end subroutine Dens_Infos

end subroutine Eval_g2_ijkl
