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
! Copyright (C) Markus P. Fuelscher                                    *
!               1992,1994, Per Ake Malmqvist                           *
!               2002, Roland Lindh                                     *
!***********************************************************************

!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
subroutine FOCKTWO(NSYM,NBAS,NFRO,KEEP,W_DLT,W_DSQ,W_FLT,nFlt,W_FSQ,LBUF,X1,X2,ExFac)

use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, DSBA_Type, Deallocate_DT
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSYM, NBAS(8), NFRO(8), KEEP(8), nFlt, LBUF
real(kind=wp), intent(in) :: W_DLT(*), W_DSQ(*), ExFac
real(kind=wp), intent(inout) :: W_FLT(nFlt), W_FSQ(*)
real(kind=wp), intent(out) :: X1(*), X2(*)
integer(kind=iwp) :: IB, IJ, IJB, IJS, IK, IOPT, IP, IPQ, IRC, IS, ISX, ISYM, JB, JK, JQ, JS, KB, KK, KLB, KS, LB, LK, LPQ, LS, &
                     LSMAX, NB, NFI, NFJ, NFK, NFL, NPQ
real(kind=wp) :: TEMP
type(DSBA_Type) :: DLT, DSQ, FLT, FSQ
real(kind=wp), external :: DDOT_

!***********************************************************************
!
! This routine has been nicked from the MOTRA package. It was
! originally written by Marcus Fuelscher, and has been slightly
! modified by P-A Malmqvist 1992-12-04.
! Further modifications by RL 2002-8-30.
! Purpose: Return FLT, which is the effective one-electron
! Hamiltonian when frozen orbitals are removed. It is the
! same as a closed-shell Fock matrix computed using the
! two-electron contribution only from frozen orbitals.
! FLT is returned as symmetry-blocked lower triangles. FSQ
! contains the same matrix, as symmetry-blocked square matrices.
! DSQ and DLT are computed in the calling routine.
! It is assumed that the SEWARD integral file was opened before
! call to this routine.
!
!***********************************************************************

call Allocate_DT(DLT,nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
call Allocate_DT(FLT,nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
call Allocate_DT(DSQ,nBas,nBas,nSym,Ref=W_DSQ)
call Allocate_DT(FSQ,nBas,nBas,nSym,Ref=W_FSQ)

do IS=1,NSYM
  IB = NBAS(IS)
  IK = KEEP(IS)
  NFI = NFRO(IS)
  do JS=1,IS
    JB = NBAS(JS)
    JK = KEEP(JS)
    NFJ = NFRO(JS)
    IJS = MUL(IS,JS)
    IJB = IB*JB
    if (IS == JS) IJB = (IB*(IB+1))/2
    do KS=1,IS
      KB = NBAS(KS)
      KK = KEEP(KS)
      NFK = NFRO(KS)
      LSMAX = KS
      if (KS == IS) LSMAX = JS
      LS = MUL(IJS,KS)
      if (LS > LSMAX) cycle
      LB = NBAS(LS)
      LK = KEEP(LS)
      NFL = NFRO(LS)
      KLB = KB*LB
      if (KS == LS) KLB = (KB*(KB+1))/2
      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      if ((IK+JK+KK+LK) /= 0) cycle
      ! NO FROZEN ORBITALS?
      if ((NFI+NFJ+NFK+NFL) == 0) cycle
      ! NO BASIS FUNCTIONS?
      if ((IJB*KLB) == 0) cycle

      if ((IS == JS) .and. (IS == KS)) then
        ! CASE 1: Integrals are of symmetry type (II/II)
        ! Coulomb and exchange terms need to be accumulated
        ! Option code 1: Begin reading at first integral.
        ! NPQ: Nr of submatrices in buffer X1.
        IOPT = 1
        LPQ = 0
        IPQ = 0
        NPQ = 0
        do IP=1,IB
          do JQ=1,IP
            IPQ = IPQ+1
            LPQ = LPQ+1
            if (IPQ > NPQ) then
              call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
              if (IRC > 1) call Error(IRC)
              ! Option code 2: Continue reading at next integral.
              IOPT = 2
              IPQ = 1
            end if
            ISX = (IPQ-1)*KLB+1
            TEMP = DDOT_(KLB,X1(ISX),1,DLT%SB(IS)%A1,1)
            FLT%SB(IS)%A1(LPQ) = FLT%SB(IS)%A1(LPQ)+TEMP
            call SQUARE(X1(ISX),X2(1),1,KB,LB)
            call DGEMV_('N',KB,LB,-Half*ExFac,X2(1),KB,DSQ%SB(IS)%A2(1:,IP),1,One,FSQ%SB(IS)%A2(1:,JQ),1)
            if (IP /= JQ) then
              call DGEMV_('N',KB,LB,-Half*ExFac,X2(1),KB,DSQ%SB(IS)%A2(1:,JQ),1,One,FSQ%SB(IS)%A2(1:,IP),1)
            end if
          end do
        end do
      else if ((IS == JS) .and. (IS /= KS)) then
        ! CASE 2: Integrals are of symmetry type (II/JJ)
        ! Coulomb terms need to be accumulated only
        IOPT = 1
        LPQ = 0
        IPQ = 0
        NPQ = 0
        do IP=1,IB
          do JQ=1,IP
            IPQ = IPQ+1
            LPQ = LPQ+1
            if (IPQ > NPQ) then
              call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
              if (IRC > 1) call Error(IRC)
              IOPT = 2
              IPQ = 1
            end if
            ISX = (IPQ-1)*KLB+1
            if (NFI /= 0) then
              TEMP = DLT%SB(IS)%A1(LPQ)
              call DAXPY_(KLB,TEMP,X1(ISX),1,FLT%SB(KS)%A1,1)
            end if
            if (NFK /= 0) then
              TEMP = DDOT_(KLB,X1(ISX),1,DLT%SB(KS)%A1,1)
              FLT%SB(IS)%A1(LPQ) = FLT%SB(IS)%A1(LPQ)+TEMP
            end if
          end do
        end do
      else if ((IS == KS) .and. (JS == LS)) then
        ! CASE 3: Integrals are of symmetry type (IJ/IJ)
        ! Exchange terms need to be accumulated only
        IOPT = 1
        LPQ = 0
        IPQ = 0
        NPQ = 0
        do IP=1,IB
          do JQ=1,JB
            IPQ = IPQ+1
            LPQ = LPQ+1
            if (IPQ > NPQ) then
              call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
              if (IRC > 1) call Error(IRC)
              IOPT = 2
              IPQ = 1
            end if
            ISX = (IPQ-1)*KLB+1
            if (NFI /= 0) then
              call DGEMV_('N',LB,KB,-Half*ExFac,X1(ISX),LB,DSQ%SB(IS)%A2(1:,IP),1,One,FSQ%SB(JS)%A2(1:,JQ),1)
            end if
            if (NFJ /= 0) then
              call DGEMV_('T',LB,KB,-Half*ExFac,X1(ISX),LB,DSQ%SB(JS)%A2(1:,JQ),1,One,FSQ%SB(IS)%A2(1:,IP),1)
            end if
          end do
        end do
      end if
    end do
  end do
end do

! Accumulate the contributions
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  do IB=1,NB
    do JB=1,IB
      IJ = IB*(IB-1)/2+JB
      FLT%SB(ISYM)%A1(IJ) = FLT%SB(ISYM)%A1(IJ)+FSQ%SB(ISYM)%A2(IB,JB)
    end do
  end do
end do

call GADSum(W_FLT,nFlt)

#ifdef _DEBUGPRINT_
write(u6,'(6X,A)') 'FROZEN FOCK MATRIX IN AO BASIS:'
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB > 0) then
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FLT%SB(ISYM)%A1,NB)
  end if
end do
#endif
call Deallocate_DT(FSQ)
call Deallocate_DT(DSQ)
call Deallocate_DT(FLT)
call Deallocate_DT(DLT)

return

contains

subroutine Error(IRC)

  integer(kind=iwp), intent(in) :: IRC

  write(u6,*) ' Error return code IRC=',IRC
  write(u6,*) ' from RDORD call, in FTWOI.'
  call Abend()

end subroutine Error

end subroutine FOCKTWO
