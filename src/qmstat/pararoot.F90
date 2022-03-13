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
!> @param[in]     Ract
!> @param[in]     BetaBol
!> @param[in]     Etot
!> @param[in,out] CalledBefore
!> @param[out]    SampleThis
!***********************************************************************

subroutine ParaRoot(Ract,BetaBol,Etot,CalledBefore,SampleThis)

implicit real*8(a-h,o-z)
external Ranf
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
#include "constants.fh"
!parameter(BoltzK=1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_)
dimension iPermutation(2,MxParT)
dimension CordstTEMP(MxCen*MxPut,3)
logical CalledBefore, WeiterBitte, Accept, SampleThis
save iTemp

BoltzK = 1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_
Dum1 = 0.0d0

! If this is first time to call on this routine.

if (.not. CalledBefore) then
  iTemp = 0
  CalledBefore = .true.
end if

! See what to do.

999 continue
if (iTemp < nTemp) then
  WeiterBitte = .true.
  iTemp = iTemp+1
  write(6,*)
  write(6,*) '    Run a new temperature ensemble...',iTemp

  ! A logical variable to make parallel tempering sampling correct.
  if (iTemp == 1) then
    SampleThis = .true.
  else
    SampleThis = .false.
  end if

else
  WeiterBitte = .false.
  iTemp = 0
  write(6,*)
  write(6,*) '    Evaluate temperature ensemble interchanges.'
end if

! If we are to run a new ensemble.

if (WeiterBitte) then
  iLuStIn = 8+nStFilT(iTemp)
  iLuStUt = 16+nStFilT(iTemp)
  write(StFilIn(6:6),'(i1.1)') nStFilT(iTemp)
  write(StFilUt(6:6),'(i1.1)') nStFilT(iTemp)

  ! Collect coordinates from proper startfile.
  call Get8(Ract,Etot)

  ! Set temperature.
  BetaBol = 1.0d0/(ParaTemps(iTemp)*BoltzK)

! If we are to attempt interchanges.

else

  do iPa=1,MxParT
    iPermutation(1,iPa) = iPa
    iPermutation(2,iPa) = iPa
  end do

  ! Construct permutations, treat nTemp == 2 as special case, the others
  ! are obtained with general algorithm.
  if (nTemp == 2) then
    iPermutation(2,1) = 2
    iPermutation(2,2) = 1
    Go To 101
  end if

  PerType = Ranf(iseed)
  if (PerType < 0.5D+0) then

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
    ! Contruct permutation for even iMac
    do iPa=2,mTemp,2
      iPermutation(2,iPa) = iPermutation(1,iPa+1)
      iPermutation(2,iPa+1) = iPermutation(1,iPa)
    end do
  end if

101 continue

  ! Now attempt interchange.

  iEnsemb = 1
2001 continue
  if (iPermutation(1,iEnsemb) == iPermutation(2,iEnsemb)) then
    iEnsemb = iEnsemb+1
    Go To 2001
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
  B1 = 1.0d0/(BoltzK*T1)
  B2 = 1.0d0/(BoltzK*T2)

  ! Make the Metropolis thing.
  BigDelta = (B2-B1)*(E2-E1)
  Expe = exp(BigDelta)
  Accept = .true.
  if (Expe < 1.0D+0) then
    Expran = ranf(iseed)
    if (Expe < Expran) Accept = .false.
  end if

  if (Accept) then
    ! Ct=C2
    iLuStIn = 8+nStFilT(iPermutation(2,iEnsemb))
    iLuStUt = 16+nStFilT(iPermutation(2,iEnsemb))
    write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
    write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(2,iEnsemb))
    call Get8(R2,E2)
    do i=1,3
      do j=1,nCent*nPart
        CordstTEMP(j,i) = Cordst(j,i)
      end do
    end do
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
    do i=1,3
      do j=1,nCent*nPart
        Cordst(j,i) = CordstTEMP(j,i)
      end do
    end do
    iLuStIn = 8+nStFilT(iPermutation(1,iEnsemb))
    iLuStUt = 16+nStFilT(iPermutation(1,iEnsemb))
    write(StFilIn(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
    write(StFilUt(6:6),'(i1.1)') nStFilT(iPermutation(1,iEnsemb))
    call Put8(R2,E2,Dum1,Dum1,Dum1)
  end if

  iEnsemb = iEnsemb+2
  if (Accept) write(6,*) '            accepted!'
  if (.not. Accept) write(6,*) '            not accepted!'
  if (iEnsemb < nTemp) Go To 2001
  write(6,*)

  ! Do some stuff before exit. The reason we go back up is that this
  ! way we will collect the right coordinates from first startfile.
  ! Observe that iTemp has been reset hence we are back at the
  ! square one again.

  Go To 999

end if

return

end subroutine ParaRoot
