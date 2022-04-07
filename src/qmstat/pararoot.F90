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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  ParaRoot
!
!> @brief
!>   Manage the parallel tempering routine
!> @author A. Ohrn
!>
!> @details
!> If our system is difficult and has small transition elements
!> in the Markov chain, we can use the parallel tempering to
!> boost sampling. This routine is the root for this; it
!> mainly handles the various configurations for the different
!> temperature ensembles; also, manages the ensemble switch.
!>
!> @param[out]    Ract
!> @param[out]    BetaBol
!> @param[out]    Etot
!> @param[in,out] CalledBefore
!> @param[out]    SampleThis
!***********************************************************************

subroutine ParaRoot(Ract,BetaBol,Etot,CalledBefore,SampleThis)

use qmstat_global, only: Cordst, iLuStIn, iLuStUt, iSeed, nCent, nPart, nStFilT, nTemp, ParaTemps, StFilIn, StFilUt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auTokJ, KBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: Ract, BetaBol, Etot
logical(kind=iwp), intent(inout) :: CalledBefore
logical(kind=iwp), intent(out) :: SampleThis
integer(kind=iwp) :: iEnsemb, iPa, iTemp = 0, mTemp
real(kind=wp) :: B1, B2, BigDelta, Dum, Dum1, E1, E2, Expe, Expran, PerType, R1, R2, T1, T2
logical(kind=iwp) :: WeiterBitte, Accept
integer(kind=iwp), allocatable :: iPermutation(:,:)
real(kind=wp), allocatable :: CordstTEMP(:,:)
real(kind=wp), parameter :: BoltzK = 1.0e-3_wp*KBoltzmann/auTokJ
real(kind=wp), external :: Random_Molcas

Dum1 = Zero

! If this is first time to call on this routine.

if (.not. CalledBefore) then
  iTemp = 0
  CalledBefore = .true.
end if

! See what to do.

call mma_allocate(CordstTEMP,3,nCent*nPart,label='CordstTEMP')
call mma_allocate(iPermutation,2,nTemp,label='iPermutation')

do
  if (iTemp < nTemp) then
    WeiterBitte = .true.
    iTemp = iTemp+1
    write(u6,*)
    write(u6,*) '    Run a new temperature ensemble...',iTemp

    ! A logical variable to make parallel tempering sampling correct.
    if (iTemp == 1) then
      SampleThis = .true.
    else
      SampleThis = .false.
    end if

  else
    WeiterBitte = .false.
    iTemp = 0
    write(u6,*)
    write(u6,*) '    Evaluate temperature ensemble interchanges.'
  end if

  if (WeiterBitte) then
    ! If we are to run a new ensemble.

    iLuStIn = 8+nStFilT(iTemp)
    iLuStUt = 16+nStFilT(iTemp)
    write(StFilIn(6:6),'(i1.1)') nStFilT(iTemp)
    write(StFilUt(6:6),'(i1.1)') nStFilT(iTemp)

    ! Collect coordinates from proper startfile.
    call Get8(Ract,Etot)

    ! Set temperature.
    BetaBol = One/(ParaTemps(iTemp)*BoltzK)
    exit

  else
    ! If we are to attempt interchanges.

    do iPa=1,nTemp
      iPermutation(:,iPa) = iPa
    end do

    ! Construct permutations, treat nTemp == 2 as special case, the others
    ! are obtained with general algorithm.
    if (nTemp == 2) then
      iPermutation(2,1) = 2
      iPermutation(2,2) = 1
    else

      PerType = Random_Molcas(iSeed)
      if (PerType < Half) then

        if (mod(nTemp,2) == 1) then
          mTemp = nTemp-1
        else
          mTemp = nTemp
        end if

        ! Construct permutation for odd iMac
        do iPa=1,mTemp,2
          iPermutation(2,iPa) = iPermutation(1,iPa+1)
          iPermutation(2,iPa+1) = iPermutation(1,iPa)
        end do

      else

        mTemp = 2*((nTemp-1)/2)
        ! Construct permutation for even iMac
        do iPa=2,mTemp,2
          iPermutation(2,iPa) = iPermutation(1,iPa+1)
          iPermutation(2,iPa+1) = iPermutation(1,iPa)
        end do
      end if

    end if

    ! Now attempt interchange.

    iEnsemb = 1
    do
      if (iPermutation(1,iEnsemb) == iPermutation(2,iEnsemb)) then
        iEnsemb = iEnsemb+1
        cycle
      end if

      ! Collect energies for the permutations.
      iLuStIn = 8+nStFilT(iPermutation(1,iEnsemb))
      iLuStUt = 16+nStFilT(iPermutation(1,iEnsemb))
      write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
      write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
      call Get8(Dum,E1)
      iLuStIn = 8+nStFilT(iPermutation(2,iEnsemb))
      iLuStUt = 16+nStFilT(iPermutation(2,iEnsemb))
      write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
      write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
      call Get8(Dum,E2)
      T1 = ParaTemps(iPermutation(1,iEnsemb))
      T2 = ParaTemps(iPermutation(2,iEnsemb))
      B1 = One/(BoltzK*T1)
      B2 = One/(BoltzK*T2)

      ! Make the Metropolis thing.
      BigDelta = (B2-B1)*(E2-E1)
      Expe = exp(BigDelta)
      Accept = .true.
      if (Expe < One) then
        Expran = Random_Molcas(iSeed)
        if (Expe < Expran) Accept = .false.
      end if

      if (Accept) then
        ! Ct=C2
        iLuStIn = 8+nStFilT(iPermutation(2,iEnsemb))
        iLuStUt = 16+nStFilT(iPermutation(2,iEnsemb))
        write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
        write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
        call Get8(R2,E2)
        CordstTEMP(:,:) = Cordst
        ! C2=C1
        iLuStIn = 8+nStFilT(iPermutation(1,iEnsemb))
        iLuStUt = 16+nStFilT(iPermutation(1,iEnsemb))
        write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
        write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
        call Get8(R1,E1)
        iLuStIn = 8+nStFilT(iPermutation(2,iEnsemb))
        iLuStUt = 16+nStFilT(iPermutation(2,iEnsemb))
        write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
        write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
        call Put8(R1,E1,Dum1,Dum1,Dum1)
        ! C1=Ct
        Cordst(:,:) = CordstTEMP
        iLuStIn = 8+nStFilT(iPermutation(1,iEnsemb))
        iLuStUt = 16+nStFilT(iPermutation(1,iEnsemb))
        write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
        write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
        call Put8(R2,E2,Dum1,Dum1,Dum1)
      end if

      iEnsemb = iEnsemb+2
      if (Accept) write(u6,*) '            accepted!'
      if (.not. Accept) write(u6,*) '            not accepted!'
      if (iEnsemb >= nTemp) exit
    end do
    write(u6,*)

    ! Do some stuff before exit. The reason we go back up is that this
    ! way we will collect the right coordinates from first startfile.
    ! Observe that iTemp has been reset hence we are back at the
    ! square one again.

  end if
end do

call mma_deallocate(CordstTEMP)
call mma_deallocate(iPermutation)

return

end subroutine ParaRoot
