************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine ChoSCF_Drv(iUHF,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,
     &                FLT,FLT_ab,nFLT,ExFac,LWFSQ,LWFSQ_ab,nOcc,nOcc_ab)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     This routine calls the original ChoSCF_Drv routine (now
C     ChoSCF_Drv_) in case of Cholesky or full DF. A new driver routine
C     is called in case of local DF (LDF).
C
      Implicit None
      Integer iUHF, nSym, nFLT, LWFSQ, LWFSQ_ab
      Integer nBas(nSym), nOcc(nSym), nOcc_ab(nSym)
      Real*8  DSQ(*), DLT(*)
      Real*8  DSQ_ab(*), DLT_ab(*)
      Real*8  FLT(*), FLT_ab(*)
      Real*8  ExFac

      Logical DoLDF

      Call DecideOnLocalDF(DoLDF)
      If (DoLDF) Then
         Call LDFSCF_Drv(iUHF,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,
     &                FLT,FLT_ab,nFLT,ExFac,LWFSQ,LWFSQ_ab,nOcc,nOcc_ab)
      Else
         Call ChoSCF_Drv_(iUHF,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,
     &                FLT,FLT_ab,nFLT,ExFac,LWFSQ,LWFSQ_ab,nOcc,nOcc_ab)
      End If

      End
      SUBROUTINE CHOSCF_DRV_(iUHF,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,
     &                FLT,FLT_ab,nFLT,ExFac,LWFSQ,LWFSQ_ab,nOcc,nOcc_ab)

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "addr.fh"
      Integer nBas(nSym), MinMem(nSym),rc
      Parameter (MaxDs = 3)
      Logical DoCoulomb(MaxDs),DoExchange(MaxDs)
      Integer Lunit(8)
      Real*8 FactC(MaxDs),FactX(MaxDs),ExFac,dmpk,dFKmat
      Integer ipDLT(MaxDs),ipDSQ(MaxDs),ipFLT(MaxDs),ipFSQ(MaxDs)
      Integer ipMSQ(MaxDs),ipNocc(MaxDs),nOcc(nSym),nOcc_ab(nSym)
      Integer nnBSF(8,8),n2BSF(8,8)
      Integer nForb(8,2),nIorb(8,2),ipMOs(2),ipKLT(2)
      Integer ALGO,NSCREEN
      Logical REORD,DECO,Cho_AUfb
      Real*8 FLT(*),FLT_ab(*)
      Real*8 DSQ(*),DSQ_ab(*),DLT(*),DLT_ab(*)
      character ww*512

      Common /CHOUNIT / Lunit
      Common /CHOSCF / REORD,DECO,dmpk,dFKmat,ALGO,NSCREEN
      Common /CHOAUF / Cho_Aufb
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV

      Integer  ip_of_Work, ip_of_iWork
      External ip_of_Work, ip_of_iWork
*
C  **************************************************
        iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
C  **************************************************

      rc=0


      do i=1,8
         Lunit(i)=-1
      end do

      call getmem('nVec_a+b','Allo','Inte',ipnVec,2*nSym)
      ipnVec_ab = ipnVec + nSym

      IF(iUHF.eq.0) THEN
         nDen = 1
         DoCoulomb(1)  = .true.
         DoExchange(1) = ExFac.ne.0.0d0 ! no SCF-exchange in pure DFT
         FactC(1)      = 1.0D0
         FactX(1)      = 1.0D0*ExFac ! ExFac used for hybrid functionals

         xFac = ExFac

         ipDLT(1) = ip_of_Work(DLT(1))
         ipDSQ(1) = ip_of_Work(DSQ(1))
         ipFLT(1) = ip_of_Work(FLT(1))
         ipFSQ(1) = LWFSQ

         If (ExFac.eq.0.0d0) Then
            CALL CHO_FOCK_DFT_RED(rc,DLT,FLT)
            If (rc.ne.0) Go To 999
            goto 997
         EndIf
*
         ipNocc(1) = ip_of_iwork(nOcc(1)) ! occup. numbers
         ipMSQ(1) = mAdCMO      ! MOs coeff as specified in addr.fh


      IF (DECO) THEN !use decomposed density

       xFac = ExFac*0.5d0

       CALL set_nnBSF(nSym,nBas,nnBSF,n2BSF)
       lVdim=0
       do i=1,nSym
          lVdim = lVdim + n2BSF(i,i)
       end do

c       koff=0
c       do i=1,nSym
c          CALL CD_TESTER(rc,Work(ipDLT(1)+koff),nBas(i),.true.)
c          CALL TRIPRT('DLT',' ',Work(ipDLT(1)+koff),nBas(i))
c          koff = koff + nnBSF(i,i)
c       end do

       CALL GETMEM('choMOs','allo','real',ipVec,lVdim)
       CALL GETMEM('ddec','allo','real',ipddec,lVdim)
       call dcopy_(lVdim,Work(ipDSQ(1)),1,Work(ipddec),1)

       ipd = ipddec
       ipV = ipVec
       Do i=1,nSym
          if(nBas(i).gt.0)then
            Ymax=0.0d0
            jaa=ipd
            do ja=1,nBas(i)
               Ymax=Max(Ymax,Work(jaa))
               jaa = jaa + nBas(i) + 1
            end do
            jaa = jaa - nBas(i) - 1
            Thr = 1.0d-8*Ymax
            CALL CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                     NumV,Thr,rc)
            If (rc.ne.0) GOTO 999
            iWork(ipnVec+i-1) = NumV
            if ( NumV .ne. nOcc(i) .and. .not.Do_SpinAV
     &                             .and. .not.Cho_Aufb ) then
            write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')
     &        'Warning! The number of occupied from the '//
     &        'decomposition of the density matrix is ',numV,
     &        ' in symm. ',i, '; Expected value = ',nOcc(i),
     &        '; Max diagonal of the density in symm. ',i,
     &        ' is equal to ',Ymax
            call WarningMessage(1,ww)
            endif
          else
            iWork(ipnVec+i-1) = 0
          endif
          ipd = ipd + n2BSF(i,i)
          ipV = ipV + n2BSF(i,i)
       End Do
       CALL GETMEM('ddec','free','real',ipddec,lVdim)

       ipNocc(1) = ipnVec ! occup. numbers

       ipMSQ(1) = ipVec       ! "Cholesky" MOs

       FactX(1) = 0.5D0*ExFac ! ExFac used for hybrid functionals

      ENDIF

      Call CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,ipNocc,
     &                ALGO,REORD,MinMem,loff1)


      if (ALGO.eq.1 .and. REORD) then

        FactX(1)=0.5d0*ExFac

      Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,
     &                FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)

            If (rc.ne.0) GOTO 999

      elseif (ALGO.eq.1 .and. .not.REORD) then

        FactX(1)=0.5d0*ExFac

        CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,
     &           FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)

            If (rc.ne.0) GOTO 999

      elseif  (ALGO.eq.2 .and. DECO) then  !use decomposed density

       FactX(1) = 0.5D0*ExFac ! vectors are scaled by construction

       if (REORD)then

          Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,
     &     FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999
       else
            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                       MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999
       endif

      elseif  (ALGO.eq.2 .and. REORD) then


      ipMSQ(1) = mAdCMO      ! MOs coeff as specified in addr.fh
      FactX(1) = 1.0D0*ExFac ! MOs coeff. are not scaled

      Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,
     &     FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999

      elseif  (ALGO.eq.2 .and. .not. REORD) then

      ipMSQ(1) = mAdCMO      ! MOs coeff as specified in addr.fh
      FactX(1) = 1.0D0*ExFac ! MOs coeff. are not scaled


            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                       MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999

      elseif (ALGO.eq.3) then

          nKB=0
          Do iSym=1,nSym
             nIorb(iSym,1) = iWork(ipNocc(1)+iSym-1)
             nKB = nKB + nIorb(iSym,1)*nBas(iSym)
          End Do
          Call GetMem('Cka','Allo','Real',ipMOs(1),nKB)

          ioff1=0
          ioff2=0
          Do iSym=1,nSym
           If (nBas(iSym)*nIorb(iSym,1).ne.0) Then
             do ikk=1,nIorb(iSym,1)
                ioff3=ioff1+nBas(iSym)*(ikk-1)
                call dcopy_(nBas(iSym),Work(ipMSQ(1)+ioff3),1,
     &                     Work(ipMOs(1)+ioff2+ikk-1),nIorb(iSym,1))
             end do
           EndIf
           ioff1=ioff1+nBas(iSym)**2
           ioff2=ioff2+nIorb(iSym,1)*nBas(iSym)
           nForb(iSym,1) = 0
          End Do

c          CALL CHO_FSCF_AO(rc,ipFLT,nForb,nOcc,ipMOs,ipDLT)  ! obsolete

          CALL CHO_FSCF(rc,nDen,ipFLT,nForb,nIorb,
     &                  ipMOs,ipDLT,xFac)


          Call GetMem('Cka','Free','Real',ipMOs(1),nKB)

          If (rc.ne.0) GOTO 999

      elseif (ALGO.eq.4) then

             ipMOs(1)=ipMSQ(1)

             Do iSym=1,nSym
                nForb(iSym,1) = 0
                nIorb(iSym,1) = iWork(ipNocc(1)+iSym-1)
             End Do

             ipKLT(1) = ipFSQ(1) ! trick to use already allocated memory

             CALL CHO_LK_SCF(rc,nDen,ipFLT,ipKLT,nForb,nIorb,
     &                         ipMOs,ipDLT,FactX(1),nSCReen,dmpk,dFKmat)

          If (rc.ne.0) GOTO 999

      else
          rc=99
          write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          CALL QUIT(rc)
      endif

      IF (DECO) CALL GETMEM('choMOs','free','real',ipVec,lVdim)

      If (ALGO.lt.3.and.ExFac.ne.0.0d0) Then

         CALL CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,
     &                  ipFLT,ipFSQ)

      EndIf
C----------------------------------------------------
 997  Continue
      Call GADSum(Work(ipFLT(1)),nFLT)


      ELSE   !  UHF calculation

         nDen = 3
C ========== Assign a truth table ==================

C --- Density(1) is Dalpha + Dbeta in a LT storage
         DoCoulomb(1)  = .true.
         DoExchange(1) = .false.
         FactC(1)      = 1.0D0

C --- Density(2) is Dalpha in a SQ storage
         DoCoulomb(2)  = .false.
         DoExchange(2) = ExFac.ne.0.0d0

C --- Density(3) is Dbeta in a SQ storage
         DoCoulomb(3)  = .false.
         DoExchange(3) = ExFac.ne.0.0d0

C --- Occupation numbers
c      call get_iarray('nIsh',nOcc,nSym)
c      call get_iarray('nIsh beta',nOcc_ab,nSym)

      ipNocc(1) = ip_of_iwork(nOcc(1)) ! dummy assignement
      ipNocc(2) = ip_of_iwork(nOcc(1)) ! occup. numbers alpha MOs
      ipNocc(3) = ip_of_iwork(nOcc_ab(1)) ! occup. numbers beta MOs

C --- MO coefficients
      ipMSQ(1)  = mAdCMO      !  dummy
      ipMSQ(2)  = mAdCMO      !  alpha MOs coeff
      ipMSQ(3)  = mAdCMO_ab   !  beta  MOs coeff

C Compute the total density Dalpha + Dbeta
      CALL DAXPY_(nFLT,1.0D0,DLT(1),1,DLT_ab(1),1)

         ipDLT(1) = ip_of_Work(DLT_ab(1)) ! total density alpha+beta LT
         ipDSQ(1) = ip_of_Work(DSQ(1))    ! dummy
         ipDSQ(2) = ip_of_Work(DSQ(1))    ! alpha density SQ
         ipDSQ(3) = ip_of_Work(DSQ_ab(1)) ! beta  density SQ

         ipFLT(1) = ip_of_Work(FLT(1))    ! Coulomb (... Falpha LT)
         ipFLT(2) = ip_of_Work(FLT_ab(1))    ! (... Fbeta LT)
         ipFSQ(2) = LWFSQ     ! alpha exchange (... Falpha SQ)
         ipFSQ(3) = LWFSQ_ab  ! beta exchange (... Fbeta SQ)

        FactX(2) = 1.0D0*ExFac ! UHF SQ-density is not scaled
        FactX(3) = 1.0D0*ExFac

         If (ExFac.eq.0.0d0) Then
            CALL CHO_FOCK_DFT_RED(rc,DLT_ab,FLT)
            If (rc.ne.0) Go To 999
            goto 998
         EndIf

      IF (DECO) THEN !use decomposed density

       CALL set_nnBSF(nSym,nBas,nnBSF,n2BSF)
       lVdim=0
       do i=1,nSym
          lVdim = lVdim + n2BSF(i,i)
       end do

       CALL GETMEM('choMOs','allo','real',ipVec,lVdim)
       CALL GETMEM('choMOs_ab','allo','real',ipVec_ab,lVdim)
       CALL GETMEM('ddec','allo','real',ipddec,lVdim)
       CALL GETMEM('ddec_ab','allo','real',ipddec_ab,lVdim)

       call dcopy_(lVdim,Work(ipDSQ(2)),1,Work(ipddec),1)
       call dcopy_(lVdim,Work(ipDSQ(3)),1,Work(ipddec_ab),1)

       ipd1 = ipddec
       ipd2 = ipddec_ab
       ipV1 = ipVec
       ipV2 = ipVec_ab
       Do i=1,nSym
          if(nBas(i).gt.0)then
              Ymax=0.0d0
              jaa=ipd1-1-nBas(i)
              do ja=1,nBas(i)
                 jaa=jaa+nBas(i)+1
                 Ymax=Max(Ymax,Work(jaa))
              end do
              Thr = 1.0d-8*Ymax
              CALL CD_InCore(Work(ipd1),nBas(i),Work(ipV1),nBas(i),
     &                       NumV1,Thr,rc)
                   If (rc.ne.0) GOTO 999
                   iWork(ipnVec+i-1) = NumV1
                   if ( NumV1 .ne. nOcc(i) .and. .not.Do_SpinAV
     &                                     .and. .not.Cho_Aufb) then
               write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')
     &           'Warning! The number of occupied from the '//
     &           'decomposition of the ALPHA dens. matrix is ',numV1,
     &           ' in symm. ',i,
     &           ';Expected value = ',nOcc(i),
     &           ';Max diagonal of the alpha density in symmetry ',
     &           i,' is equal to ',Ymax
              call WarningMessage(1,ww)
                   endif
              Ymax=0.0d0
              jaa=ipd2-1-nBas(i)
              do ja=1,nBas(i)
                 jaa=jaa+nBas(i)+1
                 Ymax=Max(Ymax,Work(jaa))
              end do
              Thr = 1.0d-8*Ymax
              CALL CD_InCore(Work(ipd2),nBas(i),Work(ipV2),nBas(i),
     &                       NumV2,Thr,rc)
                   If (rc.ne.0) GOTO 999
                   iWork(ipnVec_ab+i-1) = NumV2
                   if ( NumV2 .ne. nOcc_ab(i) .and. .not.Do_SpinAV
     &                                        .and. .not.Cho_Aufb) then
               write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')
     &           'Warning! The number of occupied from the '//
     &           'decomposition of the BETA dens. matrix is ',numV2,
     &           ' in symm. ',i,
     &           ';Expected value = ',nOcc_ab(i),
     &           ';Max diagonal of the beta density in symmetry ',
     &           i,' is equal to ',Ymax
              call WarningMessage(1,ww)
                   endif
          else
            iWork(ipnVec+i-1) = 0
            iWork(ipnVec_ab+i-1) = 0
          endif
          ipd1 = ipd1 + n2BSF(i,i)
          ipd2 = ipd2 + n2BSF(i,i)
          ipV1 = ipV1 + n2BSF(i,i)
          ipV2 = ipV2 + n2BSF(i,i)
       End Do

       CALL GETMEM('ddec_ab','free','real',ipddec_ab,lVdim)
       CALL GETMEM('ddec','free','real',ipddec,lVdim)

       ipNocc(1) = ipnVec ! dummy
       ipNocc(2) = ipnVec ! alpha occup. numbers
       ipNocc(3) = ipnVec_ab ! beta occup. numbers

       ipMSQ(1) = ipVec          ! dummy
       ipMSQ(2) = ipVec          ! "Cholesky" alpha MOs
       ipMSQ(3) = ipVec_ab       ! "Cholesky" beta  MOs

      ENDIF

      Call CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,ipNocc,
     &                ALGO,REORD,MinMem,loff1)


      if (ALGO.eq.1.and.REORD) then

      Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,
     &                FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)

            If (rc.ne.0) GOTO 999

      elseif (ALGO.eq.1 .and. .not.REORD) then

        CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,
     &           FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,ipNocc,MinMem)

            If (rc.ne.0) GOTO 999

      elseif  (ALGO.eq.2 .and. DECO) then !use decomposed density

       ipMSQ(1) = ipVec          ! dummy
       ipMSQ(2) = ipVec          ! "Cholesky" alpha MOs
       ipMSQ(3) = ipVec_ab       ! "Cholesky" beta  MOs

       if (REORD)then

          Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,
     &     FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999
       else

          CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                       MinMem,ipMSQ,ipNocc)

            If (rc.ne.0) GOTO 999
       endif

      elseif  (ALGO.eq.2 .and. REORD) then

      ipMSQ(1)  = mAdCMO      !  dummy
      ipMSQ(2)  = mAdCMO      !  alpha MOs coeff
      ipMSQ(3)  = mAdCMO_ab   !  beta  MOs coeff


      Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,
     &     FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,MinMem,ipMSQ,ipNocc)

           If (rc.ne.0) GOTO 999

      elseif  (ALGO.eq.2 .and. .not. REORD) then

      ipMSQ(1)  = mAdCMO      !  dummy
      ipMSQ(2)  = mAdCMO      !  alpha MOs coeff
      ipMSQ(3)  = mAdCMO_ab   !  beta  MOs coeff

            CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,
     &                       lOff1,FactC,FactX,ipDLT,ipDSQ,ipFLT,ipFSQ,
     &                       MinMem,ipMSQ,ipNocc)

           If (rc.ne.0) GOTO 999

      elseif (ALGO.eq.3) then

          If (DECO) Then
             ipMSQ(1) = ipVec
             ipMSQ(2) = ipVec_ab
          Else
             ipMSQ(1) = mAdCMO      ! alpha MOs coeff as in addr.fh
             ipMSQ(2) = mAdCMO_ab   ! beta MOs coeff as in addr.fh
          EndIf

          nKB1=0
          nKB2=0
          Do iSym=1,nSym
             nIorb(iSym,1) = iWork(ipNocc(2)+iSym-1)
             nIorb(iSym,2) = iWork(ipNocc(3)+iSym-1)
             nKB1 = nKB1 + nIorb(iSym,1)*nBas(iSym)
             nKB2 = nKB2 + nIorb(iSym,2)*nBas(iSym)
          End Do
          Call GetMem('Cka','Allo','Real',ipMOs(1),nKB1)
          Call GetMem('Ckb','Allo','Real',ipMOs(2),nKB2)

          ioff1=0
          ioff2=0
          ioff_ab=0
          Do iSym=1,nSym
           If (nBas(iSym)*nIorb(iSym,1).ne.0) Then
             do ikk=1,nIorb(iSym,1)
                ioff3=ioff1+nBas(iSym)*(ikk-1)
                call dcopy_(nBas(iSym),Work(ipMSQ(1)+ioff3),1,
     &                     Work(ipMOs(1)+ioff2+ikk-1),nIorb(iSym,1))
             end do
           EndIf
           If (nBas(iSym)*nIorb(iSym,2).ne.0) Then
             do ikk=1,nIorb(iSym,2)
                ioff3=ioff1+nBas(iSym)*(ikk-1)
                call dcopy_(nBas(iSym),Work(ipMSQ(2)+ioff3),1,
     &                     Work(ipMOs(2)+ioff_ab+ikk-1),nIorb(iSym,2))
             end do
           EndIf
           ioff1=ioff1+nBas(iSym)**2
           ioff2=ioff2+nIorb(iSym,1)*nBas(iSym)
           ioff_ab=ioff_ab+nIorb(iSym,2)*nBas(iSym)
           nForb(iSym,1) = 0
           nForb(iSym,2) = 0
          End Do

          nMat=2  ! alpha and beta Fock matrices

          CALL CHO_FSCF(rc,nMat,ipFLT,nForb,nIorb,
     &                  ipMOs,ipDLT,ExFac)


          Call GetMem('Cka','Free','Real',ipMOs(1),nKB1)
          Call GetMem('Ckb','Free','Real',ipMOs(2),nKB2)

          If (rc.ne.0) GOTO 999


      elseif (ALGO.eq.4) then

             nMat=2  ! alpha and beta Fock matrices

             ipMOs(1)=ipMSQ(2)
             ipMOs(2)=ipMSQ(3)

             Do iSym=1,nSym
                nForb(iSym,1) = 0
                nForb(iSym,2) = 0
                nIorb(iSym,1) = iWork(ipNocc(2)+iSym-1)
                nIorb(iSym,2) = iWork(ipNocc(3)+iSym-1)
             End Do

             ipKLT(1) = ipFSQ(2) ! trick to use already allocated memory
             ipKLT(2) = ipFSQ(3)

             CALL CHO_LK_SCF(rc,nMat,ipFLT,ipKLT,nForb,nIorb,
     &                         ipMOs,ipDLT,FactX(2),nSCReen,dmpk,dFKmat)

          If (rc.ne.0) GOTO 999

      else
          rc=99
          write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          CALL QUIT(rc)
      endif

      IF(DECO) CALL GETMEM('choMOs_ab','free','real',ipVec_ab,lVdim)
      IF(DECO) CALL GETMEM('choMOs','free','real',ipVec,lVdim)

998   Continue

C --- To get the Fbeta in LT storage ----

      If (ALGO.lt.3 .or. ExFac.eq.0.0d0) then

        call dcopy_(nFLT,Work(ipFLT(1)),1,Work(ipFLT(2)),1)

      EndIf


C --- Accumulates Coulomb and Exchange contributions
      If (ALGO.lt.3.and.ExFac.ne.0.0d0) then
         CALL CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,
     &                  ipFLT,ipFSQ)
      Endif

C----------------------------------------------------
      Call GADSum(Work(ipFLT(1)),nFLT)
      Call GADSum(Work(ipFLT(2)),nFLT)


C --- Restore the Beta-density matrix ----
C --- Copy the lower triangular of Work(ipDSQ(3))
C --- and pack the off-diagonal elements
        icount1 = 0
        icount2 = 0
        do isym=1,nsym
           koff1=ipDSQ(3) - 1 + icount1
           koff2=ipDLT(1) - 1 + icount2
           do j=1,nBas(isym)
              do k=j,nBas(isym)
               koff3 = koff1 + nBas(isym)*(j-1) + k
               koff4 = koff2 + iTri(k,j)
               if (j.eq.k) then
                Work(koff4) = Work(koff3)
               else ! packing of the matrix
                Work(koff4) = 2.0d0*Work(koff3)
               endif
              end do
           end do
           icount1 = icount1 + nBas(isym)**2
           icount2 = icount2 + nBas(isym)*(nBas(isym) + 1)/2
        end do

      ENDIF


999   If (rc.ne.0) then
         write(6,*)'CHOSCF_DRV. Non-zero return code.'
         CALL QUIT(rc)
      endif

      call getmem('nVec_a+b','Free','Inte',ipnVec,2*nSym)


      Return
      End
