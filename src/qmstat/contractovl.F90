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
!> @param[out] SintPar  The contracted basis function overlaps with extra atom--atom weights *if* this has been requested by user,
!>                      otherwise unchanged
!> @param[in]  nBaseQ   Number of AO-basis functions in QM-region
!> @param[in]  nBaseC   Like \p nBaseQ but for solvent
!> @param[in]  N        Which solvent molecule this is
!> @param[in]  nCent    How many centers the solvent molecule has
!> @param[in]  iEl      Number of elements in QM-region
!> @param[in]  nAtomsCC How many solvent atoms
!> @param[in]  iPrint   Print level
!***********************************************************************

subroutine ContractOvl(Sint,SintPar,nBaseQ,nBaseC,N,nCent,iEl,iQ_Atoms,nAtomsCC,iPrint,Inside)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "integral.fh"
#include "WrkSpc.fh"
parameter(MxSphAng=2*MxAngqNr-1)
dimension Sint(MxBas,MxBasC), SintPar(MxBas,MxBasC)
dimension Bori(3), Cori(3), Alf(MxCont), Bet(MxCont)
dimension Conkort(MxCont), Donkort(MxCont), ContrI(MxSphAng**2)
dimension Inside(MxAt,3)
logical Inside

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
      Go To 12
    end if
    do iB1=1,nBA_Q(iA1) !The basis functions on this specific atom
      do iB2=1,nBA_C(iA2)
        iQcontB = iQcontBSAV
        do iNcB1=1,nCBoA_Q(iA1,iB1) !The basis of angular type.
          iQcontB = iQcontB+1
          iCcontB = iCcontBSAV
          Bori(1) = BasOri(1,iQcontB) !Suck-out proper coord for QM.
          Bori(2) = BasOri(2,iQcontB)
          Bori(3) = BasOri(3,iQcontB)
          iqqqQ = iQang(iQcontB) !Various integers, see qfread to understand their meaning.
          nExp1 = nPrimus(iQcontB)
          nSph1 = 2*iqqqQ-1
          do i=1,nPrimus(iQcontB) !Suck-out the proper exponents for QM-region
            Alf(i) = alfa(iQcontB,i)
            Conkort(i) = cont(iQcontB,i)
          end do
          do iNcB2=1,nCBoA_C(iA2,iB2)
            iCcontB = iCcontB+1
            Cori(1) = CasOri(1,iCcontB) !Coord. of the atoms of this solvent mol.
            Cori(2) = CasOri(2,iCcontB)
            Cori(3) = CasOri(3,iCcontB)
            iqqqC = iQn(iCcontB)
            nExp2 = mPrimus(iCcontB)
            nSph2 = 2*iqqqC-1
            do j=1,mPrimus(iCcontB) !Exponents and stuff.
              Bet(j) = beta(iCcontB,j)
              Donkort(j) = dont(iCcontB,j)
            end do
            ! Now call on the routine that computes a block of primitive
            ! integrals. So if we are integrating the np-mp overlap we
            ! compute ALL primitive p-p integrals, in the first call, then
            ! they are merely contracted. This is an economical procedure
            ! for both general and ordinary contracted basis sets since all
            ! primitve overlaps are needed at some point in the contracted
            ! overlaps, the difference between general and ordinary is that
            ! in the former primitve overlaps are needed at all instances,
            ! while in the latter primitve overlaps are needed only once.
            if ((iNcB1 == 1) .and. (iNcB2 == 1)) call OverLq(Bori,Cori,Alf,Bet,iqqqQ,iqqqC,nExp1,nExp2,iPSint,Trans)
            kaunter = 0
            do i=1,nSph2 !contract
              do j=1,nSph1
                kaunter = kaunter+1
                DaNumber = 0
                do iCC=1,nExp2
                  do iCQ=1,nExp1
                    iindex = (i-1)*nSph1*nExp1+j-1+nSph1*(iCQ-1)+nSph1*nExp1*nSph2*(iCC-1)
                    DaNumber = DaNumber+Conkort(iCQ)*Donkort(iCC)*Work(iPSint+iindex)
                  end do
                end do
                ContrI(kaunter) = DaNumber
              end do
            end do
            if (iPrint >= 30) then
              write(6,*) 'Basis',iQcontB,iCcontB
              write(6,*) 'Coord.',Bori(1),Bori(2),Bori(3)
              write(6,*) 'Coord.',Cori(1),Cori(2),Cori(3)
              write(6,*) 'Alfa',(Alf(i),i=1,nPrimus(iQcontB))
              write(6,*) 'Beta',(Bet(i),i=1,mPrimus(iCcontB))
              write(6,*) 'ConQ',(Conkort(i),i=1,nPrimus(iQcontB))
              write(6,*) 'ConC',(Donkort(i),i=1,mPrimus(iCcontB))
              write(6,*) 'Angular',iqqqQ,iqqqC
              write(6,*) '#primitive',nExp1,nExp2
              write(6,*) (ContrI(k),k=1,nSph1*nSph2)
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
        ! This vector is allocated in OverLq.
        nSize = nExp1*nExp2*nSph1*nSph2
        call GetMem('AllPrims','Free','Real',iPSint,nSize)
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
12  continue
  end do
  iQcontBSAV2 = iQcontBSAV1
  iQcontBSAV = iQcontBSAV1
  iCcontBSAV = 0
  iCcontBSAV1 = 0
  iCcontBSAV2 = 0
end do
if (iPrint >= 30) then !Optional print-out.
  write(6,*)
  write(6,*) 'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE',N/nCent
  write(6,*) 'QM-AO  SOLV-AO  OVERLAP'
  do i=1,nBaseQ
    do j=1,nBaseC
      write(6,8888) i,j,Sint(i,j)
    end do
  end do
end if
8888 format(I3,'    ',I3,'       ',F12.10)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(SintPar)
  call Unused_integer(iEl)
end if

end subroutine ContractOvl
