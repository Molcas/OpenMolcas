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
! Copyright (C) 2015, Roland Lindh                                     *
!***********************************************************************

subroutine DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,nZeta_Tot,nEta_Tot,k2data1,k2data2,nAlpha,nBeta,nGamma,nDelta,ix1,iy1,iz1,ix2, &
                  iy2,iz2,ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only,NoInts,iAnga,Coor,CoorAC,mabMin,mabMax, &
                  mcdMin,mcdMax,nijkl,nabcd,mabcd,Wrk,iW2,iW4,nWork2,mWork2,HMtrxAB,HMtrxCD,la,lb,lc,ld,iCmp,iShll,NoPInts,Dij, &
                  mDij,Dkl,mDkl,Do_TnsCtl,kabcd,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl)
!***********************************************************************
! Routine for the computation of primitive integrals and accumulation  *
! to the (ab|cd) or the (e0|f0) set of integrals. If the primitive     *
! set of integrals is smaller than the set of contracted integrals     *
! the code selects to apply the HRR recursion {e0|f0} -> {ab|cd} here  *
! before the contraction generating the (ab|cd) set of integrals, if   *
! not the {e0|f0} set is contracted to the (e0|f0) set directly and    *
! HRR recursion is applied outside this routine.                       *
!                                                                      *
! For the contraction we have that either all primitive integrals      *
! can be computed in a single step, otherwise subsets of primitive     *
! integrals are computed and accumulated to the contracted set.        *
!                                                                      *
! The Wrk array is subdivided into 2 or 3 blocks depending on if the   *
! calling code iterates over subsets of primitive integrals.           *
!                                                                      *
! Memory blocking                                                      *
! ===============                                                      *
! For an iterative use:                                                *
!      iW4 points to the start of Wrk, length nWork2-mWork2            *
!      iW2 points at nWork2-mWork+1, length mWork2                     *
!      iW3 points at nWork2, length nWork3                             *
!                                                                      *
! For single iteration use:                                            *
!      iW4 and iW2 point at the start of Wrk, length nWork2            *
!      iW3 points at nWork2, length nWork3                             *
!                                                                      *
! Usage of memory                                                      *
!      Screen: does not use Wrk                                        *
!      Rys:    use iW2 section                                         *
!      HRR:    use the aggregated iW2 and iW3 section                  *
!      Cntrct: use the iW2, iW3, and iW4 sections separately           *
!                                                                      *
! Author: Roland Lindh                                                 *
!         Dept Chemistry - Angstrom, the Theoretical Chem. Prog.       *
!         Uppsala University, Uppsala, Sweden                          *
!         2015                                                         *
!***********************************************************************

use Breit, only: nComp
use k2_structure, only: k2_type
use k2_arrays, only: BraKet
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iZeta, iEta, nZeta, nEta, mZeta, mEta, nZeta_Tot, nEta_Tot, nAlpha, nBeta, nGamma, nDelta, ix1, &
                                 iy1, iz1, ix2, iy2, iz2, iAnga(4), mabMin, mabMax, mcdMin, mcdMax, nijkl, nabcd, mabcd, iW2, iW4, &
                                 nWork2, mWork2, la, lb, lc, ld, iCmp(4), iShll(4), mDij, mDkl, iBasi, jBasj, kBask, lBasl
type(k2_type), intent(in) :: k2data1, k2data2
real(kind=wp), intent(in) :: ThrInt, CutInt, vij, vkl, vik, vil, vjk, vjl, Coor(3,4), CoorAC(3,2), HMtrxAB(*), HMtrxCD(*), &
                             Dij(mDij), Dkl(mDkl), Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj), Coeff3(nGamma,kBask), &
                             Coeff4(nDelta,lBasl)
logical(kind=iwp), intent(in) :: Prescreen_On_Int_Only
logical(kind=iwp), intent(inout) :: NoInts, NoPInts, Do_TnsCtl
real(kind=wp), intent(inout) :: Wrk(nWork2)
integer(kind=iwp), intent(out) :: kabcd
integer(kind=iwp) :: i_Int, iOffE, iOffZ, iW3, lEta, lZeta, n1, n2, n3, n4, nW2, nWork3
logical(kind=iwp), parameter :: Nospecial = .false.
external :: TERI, ModU2, vCff2D, vRys2D

#ifdef _DEBUGPRINT_
write(u6,*) 'Enter DrvRys'
write(u6,*) 'iZeta, nZeta, mZeta, nZeta_Tot=',iZeta,nZeta,mZeta,nZeta_Tot
write(u6,*) 'iEta , nEta , mEta , nEta_Tot=',iEta,nEta,mEta,nEta_Tot
call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
call RecPrt('KappAB',' ',BraKet%KappaAB,1,nZeta)
call RecPrt('KappCD',' ',BraKet%KappaCD,1,nEta)
#endif

! Transfer k2 data and prescreen
! In case of integral according to Breit we still will do the prescreening according to the conventional 1/r integrals.

iOffZ = mDij-nZeta
iOffE = mDkl-nEta
call Screen(iZeta-1,iEta-1,nZeta,nEta,mZeta,mEta,lZeta,lEta, &
            k2data1,k2data2, &
            BraKet%Zeta,BraKet%ZInv,BraKet%P,BraKet%KappaAB,BraKet%IndZet, &
            BraKet%Eta,BraKet%EInv,BraKet%Q,BraKet%KappaCD,BraKet%IndEta, &
            Dij(iOffZ),Dkl(iOffE),ix1,iy1,iz1,ix2,iy2,iz2,ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only)
!write(u6,*) 'lZeta,lEta:',lZeta,lEta
if (lZeta*lEta == 0) then
  Wrk(iW2:iW2+mWork2-1) = Zero
else
  NoInts = .false.

  ! Compute [a0|c0], ijkl,a,c

  call Rys(iAnga,lZeta*lEta,BraKet%Zeta,BraKet%ZInv,lZeta,BraKet%Eta,BraKet%EInv,lEta,BraKet%P,nZeta,BraKet%Q,nEta,BraKet%KappaAB, &
           BraKet%KappaCD,Coor,Coor,CoorAC,mabMin,mabMax,mcdMin,mcdMax,Wrk(iW2),mWork2,TERI,ModU2,vCff2D,vRys2D,NoSpecial)

  ! Select between HRR before contraction or to contract
  ! and perform the HRR later once the complete set of
  ! contracted integrals have been generated.
# ifdef _DEBUGPRINT_
  write(u6,*) 'lZeta*lEta,nComp*mabcd=',lZeta*lEta,nComp,mabcd
  call RecPrt('DrvRys: [a0|c0]',' ',Wrk(iW2),lZeta*lEta,nComp*mabcd)
# endif

  if ((lZeta*lEta < nijkl) .and. (mZeta == nZeta_tot) .and. (mEta == nEta_tot)) then

    ! Apply the HRR recursions first. Note that this is only
    ! executed if used in single iteration mode. Hence,
    ! iW2 and iW4 are identical.

    n1 = lZeta*lEta*nComp*mabcd
    iW3 = iW2+n1
    call DGeTMO(Wrk(iW2),lZeta*lEta*nComp,lZeta*lEta*nComp,mabcd,Wrk(iW3),mabcd)
    Wrk(iW2:iW2+n1-1) = Wrk(iW3:iW3+n1-1)
    call TnsCtl(Wrk(iW2),nWork2,lZeta*lEta*nComp,mabMax,mabMin,mcdMax,mcdMin,HMtrxAB,HMtrxCD,la,lb,lc,ld,iCmp(1),iCmp(2), &
                iCmp(3),iCmp(4),iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
    n2 = lZeta*lEta*nComp*nabcd
    if (i_Int /= iW2) Wrk(iW2:iW2+n2-1) = Wrk(i_Int:i_Int+n2-1)
    Do_TnsCtl = .false.
    n1 = 1
    n2 = iCmp(1)*iCmp(2)
    n3 = 1
    n4 = iCmp(3)*iCmp(4)
    kabcd = nabcd
  else

    ! Postpone application of the HRR recursions until later.

    Do_TnsCtl = .true.
    n1 = mabMin
    n2 = mabMax
    n3 = mcdMin
    n4 = mcdMax
    kabcd = mabcd
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'lZeta*lEta,nComp*mabcd=',lZeta*lEta,nComp,kabcd
  call RecPrt('[a0|c0]',' ',Wrk(iW2),lZeta*lEta,nComp*kabcd)
# endif

  ! Accumulate to the contracted integrals

  if (iW4 /= iW2) then
    ! Account for size of the integrals in
    nW2 = lZeta*lEta*nComp*kabcd
  else ! iW4 == iW2
    ! Account for size of the integrals in and out
    nW2 = max(iBasi*jBasj*kBask*lBasl,lZeta*lEta)*nComp*kabcd
  end if
  iW3 = iW2+nW2
  nWork3 = mWork2-nW2
  !write(u6,*) 'iW4,iW2,iW3:',iW4,iW2,iW3
  !write(u6,*) 'nWork3:',nWork3
  call Cntrct(NoPInts,Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl,Wrk(iW2),n1,n2,n3,n4, &
              Wrk(iW3),nWork3,Wrk(iW4),BraKet%IndZet,lZeta,BraKet%IndEta,lEta,nComp)
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'nComp,kabcd,iBasi*jBasj*kBask*lBasl=',nComp,kabcd,iBasi*jBasj*kBask*lBasl
call RecPrt('DrvRys:(e0|0f)',' ',Wrk(iW4),nComp*kabcd,iBasi*jBasj*kBask*lBasl)
#endif

return

end subroutine DrvRys
