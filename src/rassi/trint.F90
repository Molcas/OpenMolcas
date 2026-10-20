!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine TRINT(CMO1,CMO2,ECORE,NGAM1,FOCKMO,NGAM2,TUVX)
!****************************************************************
! CALCULATE AND RETURN ECORE, FOCKMO, AND TUVX. ECORE IS THE
! INACTIVE-INACTIVE ENERGY CONTRIBUTION, WHICH IS TO BE MULTI-
! PLIED WITH AN OVERLAP TO GIVE A CONTRIBUTION TO THE HAMILTONIAN
! MATRIX ELEMENT. SIMILARLY, THE INACTIVE-ACTIVE CONTRIBUTION IS
! GIVEN BY THE FOCKMO ARRAY CONTRACTED WITH AN ACTIVE TRANSITION
! DENSITY MATRIX, AND THE ACTIVE-ACTIVE IS THE TRANSFORMED
! INTEGRALS IN ARRAY TUVX CONTRACTED WITH THE TRANSITION DENSITY
! TWO-ELECTRON MATRIX. THEREFORE, THE STORAGE OF THE FOCKMO AND
! THE TUVX MATRICES ARE IN THE SAME FORMAT AS THE DENSITY MATRICES.
!****************************************************************

use Index_Functions, only: nTri_Elem
use Fock_util_global, only: Fake_CMO2
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Cntrl, only: ALGO, dmpk, ERFNUC, LuOrd, Nscreen, RFPert
use Symmetry_Info, only: nIrrep
use rassi_data, only: CHFRACMEM, NASH, NASHT, NBASF, NBSQ, NBTRI, NCMO, NISH, NOSH
#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO1(NCMO), CMO2(NCMO)
integer(kind=iwp), intent(in) :: NGAM1, NGAM2
real(kind=wp), intent(out) :: ECORE, FOCKMO(NGAM1)
real(kind=wp), intent(inout) :: TUVX(NGAM2)
integer(kind=iwp) :: I, IERR, IKK, IOFF, IOFF1, IOFF2, IOFF3, iOpt, IRC, ISTA, ISTC, ISTFMO, iSym, JKK, KEEP(8), NA, nAux(8), NB, &
                     NB1, NB2, NBSX(8), nDINAO, NFAO, NI, NO, NPROD, nSymX
real(kind=wp) :: ECORE1, ECORE2
#ifdef _MOLCAS_MPP_
real(kind=wp) :: SCX
#endif
logical(kind=iwp) :: DoCholesky, FoundTwoEls, IfTest, ISQARX
type(DSBA_Type) :: Ash(2), DInAO, DLT(1), FAO, FLT(1), KSQ, MO1(2), MO2(2), Temp_SQ
real(kind=wp), allocatable :: Prod(:)
real(kind=wp), external :: DDot_

#ifdef _DEBUGPRINT_
IfTest = .true.
#else
IfTest = .false.
#endif
! THE FOLLOWING PROGRAMS USE THE ORDERED INTEGRAL FILE FOR BOTH
! THE FOCKMO MATRIX AND THE TUVX ARRAY. THEREFORE, IT IS IMPERATIVE
! THAT SYMMETRY BLOCKS HAVE NOT BEEN EXCLUDED IN THE GENERATION OF
! THIS ORDERED INTEGRAL FILE, UNLESS BOTH ACTIVE AND INACTIVE ORBITALS
! ARE MISSING FOR THAT SYMMETRY LABEL. THIS IS CHECKED FIRST:

! OPEN THE ELECTRON REPULSION INTEGRAL FILE
call f_Inquire('ORDINT',FoundTwoEls)
!call DecideOnDirect(.false.,FoundTwoEls,DoDirect,DoCholesky)
call DecideOnCholesky(DoCholesky)

if (.not. DoCholesky) then

  IOPT = 0
  IRC = 0
  call OPNORD(IRC,IOPT,'ORDINT',LUORD)
  if (IRC /= 0) then
    write(u6,*) ' ERROR: RETURN CODE=',IRC
    write(u6,*) ' RECEIVED WHEN OPENING ORDINT.'
    call ABEND()
  end if

  IRC = 0
  call GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
  if (IRC /= 0) then
    write(u6,*) ' ERROR: RETURN CODE=',IRC
    write(u6,*) ' RECEIVED WHEN CALLING GETORD.'
    call ABEND()
  end if
  if (NSYMX /= nIrrep) then
    write(u6,*)
    write(u6,*) '     *** ERROR IN SUBROUTINE TRINT ***'
    write(u6,*) '     INCOMPATIBLE NUMBERS OF IRRED. REP.'
    write(u6,*)
    call ABEND()
  end if
  do ISYM=1,nIrrep
    NB1 = NBASF(ISYM)
    NB2 = NBSX(ISYM)
    if (NB1 /= NB2) then
      write(u6,*)
      write(u6,*) '     *** ERROR IN SUBROUTINE TRINT ***'
      write(u6,*) '   INCOMPATIBLE NUMBERS OF BASIS FUNCTION'
      write(u6,*)
      call ABEND()
    end if
  end do
  do I=1,nIrrep
    if ((KEEP(I) /= 0) .and. (NOSH(I) > 0)) then
      write(u6,*) ' ERROR IN KEEP PARAMETERS ON ORDINT FILE.'
      write(u6,'(A,8I5)') ' KEEP ARRAY:',KEEP(1:nIrrep)
      write(u6,'(A,8I5)') ' NASH ARRAY:',NASH(1:nIrrep)
      write(u6,*) ' A NON-ZERO KEEP PARAMETER INDICATES THAT BASIS'
      write(u6,*) ' FUNCTIONS WITH A SPECIFIC SYMMETRY LABEL HAS'
      write(u6,*) ' BEEN SKIPPED. THIS IS POSSIBLE ONLY IF INACTIVE'
      write(u6,*) ' AND ACTIVE ORBITALS OF THAT SYMMETRY ARE MISSING.'
      call ABEND()
    end if
  end do

end if

! CALCULATE AN INACTIVE TRANSITION DENSITY MATRIX IN AO BASIS:
NDINAO = NBSQ
call Allocate_DT(DInAO,nBasF,nBasF,nIrrep)
call DIMAT(CMO1,CMO2,DINAO%A0)

NFAO = NBSQ
call Allocate_DT(FAO,nBasF,nBasF,nIrrep)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. DoCholesky) then     ! Conventional integrals
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  if (IfTest) call dVcPrt('Done',' ',DINAO%A0,NDINAO)
  ! GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
  ! PUT IT INTO A FOCK MATRIX IN AO BASIS:
  ! Note: GETH1 also adds the reaction field contribution to the
  ! 1-electron hamiltonian, and the variable ERFNuc in common /general/,
  ! which is the RF contribution to the nuclear repulsion
  call GETH1_RASSI(FAO%A0)
  if (IfTest) call dVcPrt('h0',' ',FAO%A0,NFAO)
  ! ONE CONTRIBUTION TO ECORE MUST BE CALCULATED FROM THE NAKED
  ! ONE-EL HAMILTONIAN:
  ECORE1 = DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)
  if (IfTest) write(u6,*) '      ECore1=',ECORE1

  ! ADD IN THE TWO-ELECTRON CONTRIBUTIONS TO THE FOCKAO MATRIX:

  call Allocate_DT(Temp_SQ,nBasF,nBasF,nIrrep)
  Temp_SQ%A0(:) = Zero
  call FOCK_RASSI(DINAO%A0,Temp_SQ%A0)

  ! --- FAO already contains the one-electron part
  FAO%A0(:) = FAO%A0(:)+Temp_SQ%A0(:)
  call Deallocate_DT(Temp_SQ)

# ifdef _DEBUGPRINT_
  do i=1,nIrrep
    call CHO_OUTPUT(FAO%SB(i)%A2,1,nBasF(i),1,nBasF(i),nBasF(i),nBasF(i),1,u6)
  end do
# endif

  ECORE2 = DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else       ! RI/CD integrals
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Initialize Cholesky information

  call CHO_X_INIT(irc,ChFracMem)
  if (irc /= 0) then
    write(u6,*) 'TrInt: Cho_X_Init returns error code ',irc
    call AbEnd()
  end if

  call Allocate_DT(DLT(1),nBasF,nBasF,nIrrep,aCase='TRI')
  call Fold_Mat(nIrrep,nBasF,DINAO%A0,DLT(1)%A0)

# ifdef _DEBUGPRINT_

  do i=1,nIrrep
    call cho_output(DInAO%SB(i)%A2,1,nBasF(i),1,nBasF(i),nBasF(i),nBasF(i),1,u6)
    call triprt('DLT','',DLT(1)%SB(i)%A1,nBasF(i))
  end do

# endif

  call Allocate_DT(FLT(1),nBasF,nBasF,nIrrep,aCase='TRI')

  ! GET THE ONE-ELECTRON HAMILTONIAN MATRIX FROM ONE-EL FILE AND
  ! PUT IT INTO A FOCK MATRIX IN AO BASIS:
  ! Note: CHO_GETH1 also adds the reaction field contribution to the
  ! 1-electron hamiltonian, and the variable ERFNuc in common /general/,
  ! which is the RF contribution to the nuclear repulsion

  call CHO_GETH1(nBtri,FLT(1)%A0,RFpert,ERFNuc)

  ECORE1 = DDOT_(nBtri,FLT(1)%A0,1,DLT(1)%A0,1)
  if (IfTest) write(u6,*) '      ECore1=',ECORE1,ALGO
  if (IfTest) write(u6,*) '      FAKE_CMO2=',FAKE_CMO2

# ifdef _MOLCAS_MPP_
  if (nProcs > 1) then
    scx = One/real(nProcs,kind=wp)
    ! to avoid double counting when using gadgop
    FLT(1)%A0(:) = scx*FLT(1)%A0(:)
  end if
# endif

  FAO%A0(:) = Zero ! Used as Exchange F matrix
  !                                                                    !
  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
  !                                                                    !
  if (ALGO == 1) then
    !                                                                  !
    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
    !                                                                  !
    ! reorder the MOs to fit Cholesky needs

    call Allocate_DT(MO1(1),nIsh,nBasF,nIrrep)
    call Allocate_DT(MO1(2),nIsh,nBasF,nIrrep)
    call Allocate_DT(MO2(1),nAsh,nBasF,nIrrep)
    call Allocate_DT(MO2(2),nAsh,nBasF,nIrrep)

    ioff = 0
    do iSym=1,nIrrep

      do ikk=1,nIsh(iSym)

        ioff1 = ioff+nBasF(iSym)*(ikk-1)

        MO1(1)%SB(iSym)%A2(ikk,:) = CMO1(ioff1+1:ioff1+nBasF(iSym))

        MO1(2)%SB(iSym)%A2(ikk,:) = CMO2(ioff1+1:ioff1+nBasF(iSym))
      end do

      ioff2 = ioff+nBasF(iSym)*nIsh(iSym)

      do ikk=1,nAsh(iSym)

        ioff3 = ioff2+nBasF(iSym)*(ikk-1)

        MO2(1)%SB(iSym)%A2(ikk,:) = CMO1(ioff3+1:ioff3+nBasF(iSym))

        MO2(2)%SB(iSym)%A2(ikk,:) = CMO2(ioff3+1:ioff3+nBasF(iSym))
      end do

      ioff = ioff+nBasF(iSym)*(nIsh(iSym)+nAsh(iSym))

    end do

    ! Add the two-electron contribution to the Fock matrix
    !     and compute the (tu|vx) integrals

    if (Fake_CMO2) then
      call CHO_FOCK_RASSI(DLT,MO1,MO2,FLT,TUVX,NGAM2)
    else
      call CHO_FOCK_RASSI_X(DLT,MO1,MO2,FLT,FAO,TUVX,NGAM2)
    end if

    call Deallocate_DT(MO2(2))
    call Deallocate_DT(MO2(1))
    call Deallocate_DT(MO1(2))
    call Deallocate_DT(MO1(1))
    !                                                                  !
    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
    !                                                                  !
  else  ! algo=2 (local exchange algorithm)
    !                                                                  !
    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
    !                                                                  !

    nAux(:) = nIsh(:)+nAsh(:)
    call Allocate_DT(MO1(1),nBasF,nAux,nIrrep,Ref=CMO1)
    call Allocate_DT(MO1(2),nBasF,nAux,nIrrep,Ref=CMO2)

    ! *** Only the active orbitals MO coeff need reordering
    call Allocate_DT(Ash(1),nAsh,nBasF,nIrrep)
    call Allocate_DT(Ash(2),nAsh,nBasF,nIrrep)

    do iSym=1,nIrrep

      do ikk=1,nAsh(iSym)
        jkk = nIsh(iSym)+ikk

        Ash(1)%SB(iSym)%A2(ikk,:) = MO1(1)%SB(iSym)%A2(:,jkk)

        Ash(2)%SB(iSym)%A2(ikk,:) = MO1(2)%SB(iSym)%A2(:,jkk)

      end do

    end do

    ! Add the two-electron contribution to the Fock matrix
    !     and compute the (tu|vx) integrals

    if (Fake_CMO2) then

      call CHO_LK_RASSI(DLT,MO1,FLT,FAO,TUVX,nGAM2,Ash,nScreen,dmpk)
    else

      call Allocate_DT(KSQ,nBasF,nBasF,nIrrep)
      KSQ%A0(:) = Zero

      call CHO_LK_RASSI_X(DLT,MO1,FLT,KSQ,FAO,TUVX,nGAM2,Ash,nScreen,dmpk)

      call Deallocate_DT(KSQ)
    end if

    call Deallocate_DT(Ash(2))
    call Deallocate_DT(Ash(1))
    call Deallocate_DT(MO1(2))
    call Deallocate_DT(MO1(1))

    !                                                                  !
    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
    !                                                                  !
  end if
  !                                                                    !
  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()(!
  !                                                                    !

  if (Fake_CMO2) then
    do i=1,nIrrep
      call SQUARE(FLT(1)%SB(i)%A1,FAO%SB(i)%A2,1,nBasF(i),nBasF(i))
    end do
  end if

  call Deallocate_DT(FLT(1))
  call Deallocate_DT(DLT(1))

  call GADGOP(FAO%A0,NBSQ,'+')
  call GADGOP(TUVX,NGAM2,'+')

  ECORE2 = DDOT_(NBSQ,FAO%A0,1,DINAO%A0,1)

  ! Finalize Cholesky information

  call CHO_X_FINAL(irc)
  if (irc /= 0) then
    write(u6,*) 'TrInt: Cho_X_Final returns error code ',irc
    write(u6,*) 'Try recovery -- continue.'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (IfTest) write(u6,*) '      Etwo  =',ECORE2
ECORE = Half*(ECORE1+ECORE2)
if (IfTest) write(u6,*) '      Ecore =',ECORE
! (NOTE COMPENSATION FOR DOUBLE-COUNTING OF TWO-ELECTRON CONTRIBUTION.)
call Deallocate_DT(DInAO)

! TRANSFORM THE FOCK MATRIX TO MO BASIS:
NPROD = 0
do I=1,nIrrep
  NPROD = max(NPROD,NASH(I)*NBASF(I))
end do

call mma_allocate(Prod,nProd,Label='Prod')
PROD(:) = Zero
ISTFMO = 1
ISTC = 1
FOCKMO(:) = Zero
do ISYM=1,nIrrep
  NA = NASH(ISYM)
  NI = NISH(ISYM)
  NO = NI+NA
  NB = NBASF(ISYM)
  if (NA /= 0) then
    ISTA = ISTC+NI*NB
    ! ISTFMO: BEGINNING OF F(MO,MO) BLOCK OF SYMMETRY ISYM.
    ! ISTC: BEGINNING OF MO-S OF SYMMETRY ISYM.
    ! ISTA: BEGINNING OF ACTIVE MO-S OF SYMMETRY ISYM.
    ! MATRIX MULT. PROD(AO,ACTIVE MO)=F(AO,AO)*CMO2(AO,ACTIVE MO)
    call DGEMM_('N','N',NB,NA,NB,One,FAO%SB(ISYM)%A2,NB,CMO2(ISTA),NB,Zero,PROD,NB)
    ! MATRIX MULT. F(ACT MO,ACT MO)=CMO1(AO,ACT MO)(TRP)*PROD(AO,ACT MO)
    call DGEMM_('T','N',NA,NA,NB,One,CMO1(ISTA),NB,PROD,NB,Zero,FOCKMO(ISTFMO),NASHT)
  end if
  ISTFMO = ISTFMO+NA*(NASHT+1)
  ISTC = ISTC+NO*NB
end do
call mma_deallocate(Prod)
call Deallocate_DT(FAO)

if (.not. DoCholesky) then
  ! TRANSFORM TWO-ELECTRON INTEGRALS:
  call TRAINT(CMO1,CMO2,NGAM2,TUVX)

  IRC = 0
  call CLSORD(IRC)

end if
call Chk4NaN(nTri_Elem(nasht),TUVX,iErr)
if (iErr /= 0) then
  write(u6,*) 'TrInt: TUVX corrupted'
  call Abend()
end if
call Chk4NaN(nGam1,FOCKMO,iErr)
if (iErr /= 0) then
  write(u6,*) 'TrInt: FOCKMO corrupted'
  call Abend()
end if
!call triprt('tuvx',' ',TUVX,nasht)

end subroutine TRINT
