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
      SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,W_DSQ,W_DLT,W_FLT,ExFac,FSQ, &
     &                      W_CMO)

      use Data_Structures, only: DSBA_Type, Allocate_DSBA,              &
     &                           Deallocate_DSBA, Integer_Pointer
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"

      Type (DSBA_Type) FSQ

      Integer nBas(8), MinMem(8),rc
      Real*8 W_FLT(*),W_CMO(*), W_DSQ(*),W_DLT(*)
      Parameter (MaxDs = 1)
      Logical DoCoulomb(MaxDs),DoExchange(MaxDs)
      Real*8 FactC(MaxDs),FactX(MaxDs),ExFac
      Integer, Target :: nOcc(nSym)

      Integer, Allocatable, Target :: nVec(:)

      Type (Integer_Pointer) :: pNocc(1)

      Type (DSBA_Type) Vec, DDec, DLT, FLT, DSQ, CMO, MSQ(MaxDs)

#include "chounit.fh"
#include "choras.fh"
!                                                                      *
!***********************************************************************
!                                                                      *
      rc=0

      Lunit(:)=-1

      nDen = 1
      DoCoulomb(1)  = .true.
      DoExchange(1) = ExFac.ne.Zero
      FactC(1)      = One
      FactX(1)      = Half*ExFac ! ExFac used for hybrid functionals

      Call Allocate_DSBA(CMO,nBas,nBas,nSym,Ref=W_CMO)

      Call Allocate_DSBA(DLT,nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
      Call Allocate_DSBA(FLT,nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)

      Call Allocate_DSBA(DSQ,nBas,nBas,nSym,Ref=W_DSQ)

      iUHF=0

       IF (DECO) THEN !use decomposed density
! ==============  Alternative A: Use decomposed density matrix =====
       Call mma_allocate(nVec,nSym,Label='nVec')

! Allocate vectors representing decomposed density matrix:
       Call Allocate_DSBA(Vec,nBas,nBas,nSym)

! ------------------------------------------------------------------
       Call Allocate_DSBA(Ddec,nBas,nBas,nSym)
       DDec%A0(:) = DSQ%A0(:)
       Do i=1,nSym
! Loop over symmetries
          if(nBas(i).gt.0)then
            Ymax=0.0d0
            do ja=1,nBas(i)
               Ymax=Max(Ymax,DDec%SB(i)%A2(ja,ja))
            end do
            Thr = 1.0d-13*Ymax
! Call for decomposition:
            CALL CD_InCore(DDec%SB(i)%A2,nBas(i),Vec%SB(i)%A2,nBas(i),  &
     &                     NumV,Thr,rc)
            If (rc.ne.0) GOTO 999
            nVec(i) = NumV

            if ( NumV .ne. nOcc(i) ) then
               write(6,*)'Warning! The number of occupied from the dec',&
     &'omposition of the density matrix is ',numV,' in symm. ',i
               write(6,*)'Expected value = ',nOcc(i)
               write(6,*)'Max diagonal of the density in symm. ',i,     &
     &' is equal to ',Ymax
            endif

          else
            nVec(i) = 0
          endif
! End of loop over symmetries
       End Do
       Call Deallocate_DSBA(DDec)
! ------------------------------------------------------------------


       pNocc(1)%I1(1:)=>nVec(1:)

       Call Allocate_DSBA(MSQ(1),nBas,nBas,nSym,Ref=Vec%A0)

! ========End of  Alternative A: Use decomposed density matrix =====
      ELSE

       pNocc(1)%I1(1:)=>nOcc(1:)

       Call Allocate_DSBA(MSQ(1),nBas,nBas,nSym,Ref=CMO%A0)

      ENDIF

      Call CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,pNocc,                  &
     &                ALGO,REORD,MinMem,loff1)

! Here follows a long if nest with six combinations:
! ALGO is 1 ,  REORD is .true. or .false., or
! ALGO is 2,  REORD is .true. or .false., DECO is .true.or .false.
      if (ALGO.eq.1) then
        if (REORD) then
! ALGO.eq.1.and.REORD:
      Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,    &
     &                FactX,[DLT],[DSQ],[FLT],[FSQ],pNocc,MinMem)
           If (rc.ne.0) GOTO 999

        else
! ALGO.eq.1.and. .not.REORD:
        CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,         &
     &           FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],pNocc,MinMem)
           If (rc.ne.0) GOTO 999
        end if

      elseif (ALGO.eq.2) then
        if (DECO) then !use decomposed density

          FactX(1) = 0.5D0*ExFac ! vectors are scaled by construction
          if (REORD)then
! ALGO.eq.2.and.DECO.and.REORD:
            Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,    &
     &                  lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],      &
     &                  MinMem,MSQ,pNocc)
            If (rc.ne.0) GOTO 999

          else
! ALGO.eq.2.and.DECO.and. .not.REORD:
            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,              &
     &                  lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],      &
     &                  MinMem,MSQ,pNocc)
            If (rc.ne.0) GOTO 999
          endif


        else
          if (REORD) then
! ALGO.eq.2.and. ..not.DECO.and.REORD:
            FactX(1) = 1.0D0*ExFac ! because MOs coeff. are not scaled
            Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,    &
     &                lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],        &
     &                MinMem,MSQ,pNocc)
            if (rc.ne.0) GOTO 999
          else
! ALGO.eq.2.and. ..not.DECO.and.REORD:
            FactX(1) = 1.0D0*ExFac ! because MOs coeff. are not scaled
            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,              &
     &                lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],        &
     &                MinMem,MSQ,pNocc)
            if (rc.ne.0) GOTO 999
          end if
        end if

      else
        rc=99
        write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
        CALL QUIT(rc)
      endif

      CALL CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,[FLT],[FSQ])

 999    continue
      if (rc.ne.0) then
         write(6,*)'CHORAS_DRV. Non-zero return code. rc= ',rc
         CALL QUIT(rc)
      end if

      pNocc(1)%I1=>Null()
      Call Deallocate_DSBA(MSQ(1))
      IF (DECO) Then
       Call Deallocate_DSBA(Vec)
       Call mma_deallocate(nVec)
      End If
      Call Deallocate_DSBA(DSQ)
      Call Deallocate_DSBA(DLT)
      Call Deallocate_DSBA(FLT)
      Call Deallocate_DSBA(CMO)

      Return
      End
