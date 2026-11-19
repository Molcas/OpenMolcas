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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

subroutine ROMP2_Driver(ireturn)
use Definitions, only: wp,iwp,u6
use MBPT2_Global, only: MBPT2_All_Finite
use Constants, only: Zero
use UMP2_Global
use Para_Info, only: MyRank,Is_Real_Par,nProcs
use ROMP2_Reference, only: ROMP2_Read,ROMP2_Fock
use ROMP2_Cholesky, only: ROMP2_Cho_Open,ROMP2_Cho_Close,ROMP2_Cho_Energy
use ROMP2_Semicanonical, only: ROMP2_Semi
use ROMP2_Conventional, only: ROMP2_Conventional_Energy
use stdalloc, only: mma_allocate,mma_deallocate
implicit none

interface
  subroutine RdInp(CMO,Eall,Eocc,Eext,iTst,ESCF,OpenErr,OpenMessage,SymFrozen,SymOcc)
    use Definitions, only: wp,iwp
#include "intent.fh"
    real(kind=wp),optional,intent(_OUT_) :: CMO(*),Eall(*),Eocc(*),Eext(*)
    integer(kind=iwp),optional,intent(out) :: iTst,OpenErr
    real(kind=wp),optional,intent(out) :: ESCF
    character(len=*),optional,intent(out) :: OpenMessage
    integer(kind=iwp),optional,intent(out) :: SymFrozen(:)
    integer(kind=iwp),optional,intent(in) :: SymOcc(:,:)
  end subroutine RdInp
end interface

#include "warnings.h"
integer(kind=iwp),intent(out) :: ireturn
integer(kind=iwp) :: ierr,cleanupRC
real(kind=wp) :: E1A,E1B,Singles(2)
real(kind=wp),allocatable :: FA(:,:),FB(:,:)
character(len=256) :: Message
logical :: NonfiniteTotal
ireturn=_RC_INPUT_ERROR_
call UMP2_Clean()
call RunCalculation()
call ROMP2_Cho_Close(cleanupRC)
if ((cleanupRC/=0) .and. ((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) &
  write(u6,'(A,I8)') ' ROMP2 Cholesky cleanup return code: ',cleanupRC
if ((ireturn/=_RC_ALL_IS_WELL_).and.((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) write(u6,'(A)') ' ROMP2 error: '//trim(Message)
call mma_deallocate(FA,safe='*')
call mma_deallocate(FB,safe='*')
call UMP2_Clean()

contains

subroutine RunCalculation()
integer(kind=iwp) :: MPIFlag,i
call ROMP2_Read(ierr,Message)
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
call RdInp(OpenErr=ierr,OpenMessage=Message)
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
MPIFlag=0
if ((nFro<0).or.(nFro>min(nOccA,nOccB))) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  Message='ROMP2 Frozen must be between zero and the number of doubly occupied reference orbitals.'
  return
end if
if (DoCholesky) then
  call ROMP2_Cho_Open(ierr,Message)
  MPIFlag=0
  if (ierr/=0) MPIFlag=1
  if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
  if (MPIFlag/=0) then
    if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
    ierr=1
    Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
  end if
  if (ierr/=0) return
end if
call mma_allocate(FA,nBas,nBas,label='ROMP2 Fock alpha')
call mma_allocate(FB,nBas,nBas,label='ROMP2 Fock beta')
call ROMP2_Fock(FA,FB,ierr,Message)
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
! Root diagonalizes once: every rank must use the same orbital gauge.
ierr=0
E1A=Zero
E1B=Zero
if (((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0))) then
  call ROMP2_Semi(CAlpha,FA,nBas,nOccA,nFro,EOrbA,E1A,ierr)
  if (ierr/=0) then
    Message='ROMP2 alpha semicanonicalization or singles evaluation failed.'
  else
    call ROMP2_Semi(CBeta,FB,nBas,nOccB,nFro,EOrbB,E1B,ierr)
    if (ierr/=0) Message='ROMP2 beta semicanonicalization or singles evaluation failed.'
  end if
end if
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
if (((nProcs>1).and.Is_Real_Par())) then
  if (MyRank/=0) then
    CAlpha=Zero
    CBeta=Zero
    EOrbA=Zero
    EOrbB=Zero
  end if
  call GADGOp(CAlpha(1,1),size(CAlpha),'+')
  call GADGOp(CBeta(1,1),size(CBeta),'+')
  call GADGOp(EOrbA(1),size(EOrbA),'+')
  call GADGOp(EOrbB(1),size(EOrbB),'+')
  Singles(:)=[E1A,E1B]
  call GADGOp(Singles(1),2,'+')
  E1A=Singles(1)
  E1B=Singles(2)
end if
call mma_deallocate(FA)
call mma_deallocate(FB)
call UMP2_Setup_Spaces(ierr,Message)
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
if (DoCholesky) then
  call ROMP2_Cho_Energy(ierr,Message)
else
  call ROMP2_Conventional_Energy(ierr,Message)
end if
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) return
ECorr=E1A+E1B+EAA+EBB+EAB
ETotal=ESCF+ECorr
! Assignment converts the intrinsic result to the default logical kind
! expected by MBPT2_Any, including builds with promoted default kinds.
NonfiniteTotal=.not.MBPT2_All_Finite(ETotal)
MPIFlag=0
if (NonfiniteTotal) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  Message='Nonfinite ROMP2 total energy.'
  return
end if
! Finalize before publishing any energy as successful.
call ROMP2_Cho_Close(ierr)
MPIFlag=0
if (ierr/=0) MPIFlag=1
if ((nProcs>1).and.Is_Real_Par()) call gaIgOP_SCAL(MPIFlag,'max')
if (MPIFlag/=0) then
  if (ierr/=0) write(u6,'(A,I6,2A)') ' MBPT2 rank ',MyRank,': ',trim(Message)
  ierr=1
  Message='MBPT2 failed on at least one rank: see preceding diagnostics.'
end if
if (ierr/=0) then
  Message='ROMP2 Cholesky finalization failed.'
  return
end if
if (((nProcs<=1).or.(.not.Is_Real_Par()).or.(MyRank==0)) .and. (iPL >= 2)) then
  call PrintResultsHeader()

  write(u6,'(6X,A,T48,A)') 'Method', '= Semicanonical ROHF-MBPT(2)'
  if (DoCholesky) then
    write(u6,'(6X,A,T48,A)') 'Integral treatment', '= Cholesky'
  else
    write(u6,'(6X,A,T48,A)') 'Integral treatment', '= Conventional stored AO integrals'
  end if
  write(u6,'(6X,A,T48,A,I8)') 'Frozen occupied orbitals per spin','=',nFro

  if (iPL >= 3) then
    write(u6,'(6X,A,T48,A,2I8)') 'Occupied orbitals (alpha,beta)','=',nOccA,nOccB
    write(u6,'(6X,A,T48,A,2I8)') 'Correlated occupied orbitals (alpha,beta)','=',nOA,nOB
    write(u6,'(6X,A,T48,A,2I8)') 'Virtual orbitals (alpha,beta)','=',nVA,nVB

    if (nFro > 0) then
      write(u6,*)
      write(u6,'(6X,A)') &
        'Frozen-core convention: lowest occupied orbitals after spin semicanonicalization'
      write(u6,'(6X,A)') &
        'Index          alpha orbital energy / a.u.      beta orbital energy / a.u.'
      do i=1,nFro
        write(u6,'(6X,I5,2F31.12)') i,EOrbA(i),EOrbB(i)
      end do
    end if
  end if

  write(u6,'(//6X,A,T48,A,F20.10,A)') 'ROHF reference energy','=',ESCF,' a.u.'
  write(u6,'(6X,A,T48,A,F20.10,A)') 'Second-order correlation energy','=',ECorr,' a.u.'
  write(u6,'(6X,A,T48,A,F20.10,A)') 'Total ROMP2 energy','=',ETotal,' a.u.'
  write(u6,*)
  call PrintResult(u6,'(6X,A,T50,F19.10)', &
                   'Total MBPT2 energy',0,'',[ETotal],1)
  write(u6,*)
end if
EnergyReady=.true.
call Store_Energies(1,[ETotal],1)
call Put_cArray('Relax Method','ROMP2   ',8)
call Put_iScalar('mp2prpt',0)
call Add_Info('E_ROMP2',[ETotal],1,8)
call Add_Info('E_ROMP2_SINGLES',[E1A+E1B],1,8)
call Add_Info('E_ROMP2_AA',[EAA],1,8)
call Add_Info('E_ROMP2_BB',[EBB],1,8)
call Add_Info('E_ROMP2_AB',[EAB],1,8)
ireturn=_RC_ALL_IS_WELL_
end subroutine RunCalculation

subroutine PrintResultsHeader()

  write(u6,'(/,6X,A)') &
    '********************************************************************************************************'
  write(u6,'(6X,A)') &
    '*                                               Results                                                *'
  write(u6,'(6X,A/)') &
    '********************************************************************************************************'

end subroutine PrintResultsHeader

end subroutine ROMP2_Driver
