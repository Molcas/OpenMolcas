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
!  ContractOvl
!
!> @brief
!>   Compute the overlaps between solvent and solute in contracted basis-functions
!> @author A. Ohrn
!>
!> @details
!> Here the overlap between the QM-region contracted AO-basis
!> functions and the present solvent molecule contracted AO-basis
!> functions are computed. In order to use the fact that we use
!> contracted functions to the maximum, we compute the overlaps with
!> primitive functions only once, then we transform this matrix to
!> all relevant contracted overlaps. After that, the old primitive
!> integrals are discarded and a new set of primitive are computed.
!> This is very nice since ::OverLq is rather slow. The problems we
!> get are that we must use rather elaborate schemes to get right
!> digit in right place.
!>
!> @param[out] Sint     The contracted basis function overlaps
!> @param[in]  nBaseQ   Number of AO-basis functions in QM-region
!> @param[in]  nBaseC   Like \p nBaseQ but for solvent
!> @param[in]  N        Which solvent molecule this is
!> @param[in]  nCent    How many centers the solvent molecule has
!> @param[in]  iQ_Atoms
!> @param[in]  nAtomsCC How many solvent atoms
!> @param[in]  iPrint   Print level
!> @param[in]  Inside
!***********************************************************************

subroutine ContractOvl(Sint,nBaseQ,nBaseC,N,nCent,iQ_Atoms,nAtomsCC,iPrint,Inside)

use qmstat_global, only: Alfa, BasOri, Beta, CasOri, Cont, Dont, iQang, iQn, iWoGehenC, iWoGehenQ, mPrimus, MxAngqNr, nBA_C, &
                         nBA_Q, nBonA_C, nBonA_Q, nCBoA_C, nCBoA_Q, nPrimus
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nBaseQ, nBaseC, N, nCent, iQ_Atoms, nAtomsCC, iPrint
real(kind=wp), intent(_OUT_) :: Sint(nBaseQ,nBaseC)
logical(kind=iwp), intent(in) :: Inside(iQ_Atoms,nAtomsCC)
integer(kind=iwp) :: i, iA1, iA2, iB1, iB2, iC, iCC, iCcontB, iCcontBSAV, iCcontBSAV1, iCcontBSAV2, iCQ, iNcB1, iNcB2, iQ, &
                     iQcontB, iQcontBSAV, iQcontBSAV1, iQcontBSAV2, iqqqC, iqqqQ, j, kaunter, kreichner, nExp1, nExp2, nSph1, nSph2
real(kind=wp) :: Bori(3), Cori(3), DaNumber
real(kind=wp), allocatable :: Alf(:), Bet(:), Conkort(:), ContrI(:), Donkort(:), PSint(:,:,:,:)

call mma_allocate(Alf,size(Alfa,2),label='Alf')
call mma_allocate(Bet,size(Beta,2),label='Bet')
call mma_allocate(Conkort,size(Cont,2),label='Conkort')
call mma_allocate(Donkort,size(Dont,2),label='Donkort')
call mma_allocate(ContrI,(2*MxAngqNr-1)**2,label='ContrI')

nSph1 = 0
nSph2 = 0
iQcontBSAV = 0
iCcontBSAV = 0
iQcontBSAV1 = 0
iCcontBSAV1 = 0
iQcontBSAV2 = 0
iCcontBSAV2 = 0
do iA1=1,iQ_Atoms !The atoms
  do iA2=1,nAtomsCC
    if (.not. Inside(iA1,iA2)) then !when atom-pair too far from each other, do this then skip.
      iCcontBSAV = iCcontBSAV+nBonA_C(iA2)
      iCcontBSAV1 = iCcontBSAV
      iCcontBSAV2 = iCcontBSAV
      if (iA2 == nAtomsCC) then
        iQcontBSAV = iQcontBSAV+nBonA_Q(iA1)
        iQcontBSAV1 = iQcontBSAV
      end if
      cycle
    end if
    do iB1=1,nBA_Q(iA1) !The basis functions on this specific atom
      do iB2=1,nBA_C(iA2)
        iQcontB = iQcontBSAV
        do iNcB1=1,nCBoA_Q(iA1,iB1) !The basis of angular type.
          iQcontB = iQcontB+1
          iCcontB = iCcontBSAV
          Bori(:) = BasOri(:,iQcontB) !Suck-out proper coord for QM.
          iqqqQ = iQang(iQcontB) !Various integers, see qfread to understand their meaning.
          nExp1 = nPrimus(iQcontB)
          nSph1 = 2*iqqqQ-1
          ! Suck-out the proper exponents for QM-region
          Alf(1:nPrimus(iQcontB)) = alfa(iQcontB,1:nPrimus(iQcontB))
          Conkort(1:nPrimus(iQcontB)) = cont(iQcontB,1:nPrimus(iQcontB))
          do iNcB2=1,nCBoA_C(iA2,iB2)
            iCcontB = iCcontB+1
            Cori(:) = CasOri(:,iCcontB) !Coord. of the atoms of this solvent mol.
            iqqqC = iQn(iCcontB)
            nExp2 = mPrimus(iCcontB)
            nSph2 = 2*iqqqC-1
            ! Exponents and stuff.
            Bet(1:mPrimus(iCcontB)) = beta(iCcontB,1:mPrimus(iCcontB))
            Donkort(1:mPrimus(iCcontB)) = dont(iCcontB,1:mPrimus(iCcontB))
            ! Now call on the routine that computes a block of primitive
            ! integrals. So if we are integrating the np-mp overlap we
            ! compute ALL primitive p-p integrals, in the first call, then
            ! they are merely contracted. This is an economical procedure
            ! for both general and ordinary contracted basis sets since all
            ! primitive overlaps are needed at some point in the contracted
            ! overlaps, the difference between general and ordinary is that
            ! in the former primitive overlaps are needed at all instances,
            ! while in the latter primitive overlaps are needed only once.
            if ((iNcB1 == 1) .and. (iNcB2 == 1)) then
              call mma_allocate(PSint,nSph1,nExp1,nSph2,nExp2,label='AllPrims')
              call OverLq(Bori,Cori,Alf,Bet,iqqqQ,iqqqC,nExp1,nExp2,PSint)
            end if
            kaunter = 0
            do i=1,nSph2 !contract
              do j=1,nSph1
                kaunter = kaunter+1
                DaNumber = Zero
                do iCC=1,nExp2
                  do iCQ=1,nExp1
                    DaNumber = DaNumber+Conkort(iCQ)*Donkort(iCC)*PSint(j,iCQ,i,iCC)
                  end do
                end do
                ContrI(kaunter) = DaNumber
              end do
            end do
            if (iPrint >= 30) then
              write(u6,*) 'Basis',iQcontB,iCcontB
              write(u6,*) 'Coord.',Bori(1),Bori(2),Bori(3)
              write(u6,*) 'Coord.',Cori(1),Cori(2),Cori(3)
              write(u6,*) 'Alfa',(Alf(i),i=1,nPrimus(iQcontB))
              write(u6,*) 'Beta',(Bet(i),i=1,mPrimus(iCcontB))
              write(u6,*) 'ConQ',(Conkort(i),i=1,nPrimus(iQcontB))
              write(u6,*) 'ConC',(Donkort(i),i=1,mPrimus(iCcontB))
              write(u6,*) 'Angular',iqqqQ,iqqqC
              write(u6,*) '#primitive',nExp1,nExp2
              write(u6,*) ContrI(1:nSph1*nSph2)
            end if
            kreichner = 0
            do iC=1,nSph2
              do iQ=1,nSph1
                kreichner = kreichner+1
                Sint(iWoGehenQ(iQcontB,iQ),iWoGehenC(iCcontB,iC)) = ContrI(kreichner)
              end do
            end do
          end do
        end do
        iCcontBSAV = iCcontB
        call mma_deallocate(PSint)
      end do
      ! OH NO!, these things have to do with
      ! getting the right number in right
      ! place. This should probably be
      ! changed sometime to a less cumbersome method.
      iQcontBSAV = iQcontB
      iCcontBSAV1 = iCcontBSAV
      iCcontBSAV = iCcontBSAV2
    end do
    iQcontBSAV1 = iQcontBSAV
    iQcontBSAV = iQcontBSAV2
    iCcontBSAV2 = iCcontBSAV1
    iCcontBSAV = iCcontBSAV1
  end do
  iQcontBSAV2 = iQcontBSAV1
  iQcontBSAV = iQcontBSAV1
  iCcontBSAV = 0
  iCcontBSAV1 = 0
  iCcontBSAV2 = 0
end do
call mma_deallocate(Alf)
call mma_deallocate(Bet)
call mma_deallocate(Conkort)
call mma_deallocate(Donkort)
call mma_deallocate(ContrI)
if (iPrint >= 30) then !Optional print-out.
  write(u6,*)
  write(u6,*) 'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE',N/nCent
  write(u6,*) 'QM-AO  SOLV-AO  OVERLAP'
  do i=1,nBaseQ
    do j=1,nBaseC
      write(u6,8888) i,j,Sint(i,j)
    end do
  end do
end if

return

8888 format(I3,'    ',I3,'       ',F12.10)

end subroutine ContractOvl
