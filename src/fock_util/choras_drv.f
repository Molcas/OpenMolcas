************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,DSQ,DLT,FLT,
     &                      ExFac,LWFSQ,CMO)

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Integer nBas(8), MinMem(8),rc
      Real*8 FLT(*),CMO(*)
      Real*8 DSQ(*),DLT(*)
      Parameter (MaxDs = 1)
      Logical DoCoulomb(MaxDs),DoExchange(MaxDs)
      Integer Lunit(8)
      Real*8 FactC(MaxDs),FactX(MaxDs),ExFac
      Integer ipDLT(MaxDs),ipDSQ(MaxDs),ipFLT(MaxDs),ipFSQ(MaxDs)
      Integer ipMSQ(MaxDs),ipNocc(MaxDs),nOcc(nSym)
      Integer nnBSF(8,8),n2BSF(8,8)
      Integer ALGO
      Logical REORD,DECO

      Common /CHOUNIT / Lunit
      Common /CHORAS / REORD,DECO,ALGO
*
*
C  **************************************************

      rc=0

      do i=1,8
         Lunit(i)=-1
      end do

      nDen = 1
      DoCoulomb(1)  = .true.
      DoExchange(1) = ExFac.ne.0.0d0
      FactC(1)      = 1.0D0
      FactX(1)      = 0.5D0*ExFac ! ExFac used for hybrid functionals


      ipDLT(1) = ip_of_Work(DLT(1))
      ipDSQ(1) = ip_of_Work(DSQ(1))
      ipFLT(1) = ip_of_Work(FLT(1))
      ipFSQ(1) = LWFSQ

      iUHF=0

      ipNocc(1) = ip_of_iwork(nOcc(1)) ! occup. numbers

       IF (DECO) THEN !use decomposed density
* ==============  Alternative A: Use decomposed density matrix =====
       call getmem('nVec','Allo','Inte',ipnVec,nSym)

       CALL set_nnBSF(nSym,nBas,nnBSF,n2BSF)
       lVdim=0
       do i=1,nSym
          lVdim = lVdim + n2BSF(i,i)
       end do

* Allocate vectors representing decomposed density matrix:
       CALL GETMEM('choMOs','allo','real',ipVec,lVdim)

* ------------------------------------------------------------------
       CALL GETMEM('ddec','allo','real',ipddec,lVdim)
       call dcopy_(lVdim,Work(ipDSQ(1)),1,Work(ipddec),1)
       ipd = ipddec
       ipV = ipVec
       Do i=1,nSym
* Loop over symmetries
          if(nBas(i).gt.0)then
            Ymax=0.0d0
            do ja=1,nBas(i)
               jaa=ipd-1+nBas(i)*(ja-1)+ja
               Ymax=Max(Ymax,Work(jaa))
            end do
            Thr = 1.0d-13*Ymax
* Call for decomposition:
            CALL CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                     NumV,Thr,rc)
            If (rc.ne.0) GOTO 999
            iwork(ipnVec+i-1) = NumV

            if ( NumV .ne. nOcc(i) ) then
               write(6,*)'Warning! The number of occupied from the dec',
     &'omposition of the density matrix is ',numV,' in symm. ',i
               write(6,*)'Expected value = ',nOcc(i)
               write(6,*)'Max diagonal of the density in symm. ',i,
     &' is equal to ',Ymax
            endif

          else
            iwork(ipnVec+i-1) = 0
          endif
          ipd = ipd + n2BSF(i,i)
          ipV = ipV + n2BSF(i,i)
* End of loop over symmetries
       End Do
       CALL GETMEM('ddec','free','real',ipddec,lVdim)
* ------------------------------------------------------------------


       ipNocc(1) = ipnVec ! occup. numbers

       ipMSQ(1) = ipVec       ! "Cholesky" MOs

* ========End of  Alternative A: Use decomposed density matrix =====
      ENDIF

      Call CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,ipNocc,
     &                ALGO,REORD,MinMem,loff1)

* Here follows a long if nest with six combinations:
* ALGO is 1 ,  REORD is .true. or .false., or
* ALGO is 2,  REORD is .true. or .false., DECO is .true.or .false.
      if (ALGO.eq.1) then
        if (REORD) then
* ALGO.eq.1.and.REORD:
      Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,
     &                FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)
           If (rc.ne.0) GOTO 999

        else
* ALGO.eq.1.and. .not.REORD:
        CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,
     &           FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)
           If (rc.ne.0) GOTO 999
        end if

      elseif (ALGO.eq.2) then
        if (DECO) then !use decomposed density

          FactX(1) = 0.5D0*ExFac ! vectors are scaled by construction
          if (REORD)then
* ALGO.eq.2.and.DECO.and.REORD:
            Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,
     &                  lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                  MinMem,ipMSQ,ipNocc)
            If (rc.ne.0) GOTO 999

          else
* ALGO.eq.2.and.DECO.and. .not.REORD:
            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                  lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                  MinMem,ipMSQ,ipNocc)
            If (rc.ne.0) GOTO 999
          endif


        else
          if (REORD) then
* ALGO.eq.2.and. ..not.DECO.and.REORD:
            ipMSQ(1) = ip_of_work(CMO(1))
            FactX(1) = 1.0D0*ExFac ! because MOs coeff. are not scaled
            Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,
     &                lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                MinMem,ipMSQ,ipNocc)
            if (rc.ne.0) GOTO 999
          else
* ALGO.eq.2.and. ..not.DECO.and.REORD:
            ipMSQ(1) = ip_of_work(CMO(1))
            FactX(1) = 1.0D0*ExFac ! because MOs coeff. are not scaled
            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                MinMem,ipMSQ,ipNocc)
            if (rc.ne.0) GOTO 999
          end if
        end if

      else
        rc=99
        write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
        CALL QUIT(rc)
      endif

      CALL CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,ipFLT,ipFSQ)

 999    continue
      if (rc.ne.0) then
         write(6,*)'CHORAS_DRV. Non-zero return code. rc= ',rc
         CALL QUIT(rc)
      end if

      IF (DECO) Then
       CALL GETMEM('choMOs','free','real',ipVec,lVdim)
       call getmem('nVec','Free','Inte',ipnVec,nSym)
      End If

      Return
      End
