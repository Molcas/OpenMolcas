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
      Subroutine NatOrb_RASSCF(CMOO,SCR1,SCR2,SMAT,CMON,OCCN)
!
!     RASSCF program: version IBM-3090: Output section
!
!     PURPOSE: Calculation of natural orbitals from the
!              density matrix. These orbitals are written onto JOBIPH
!              in position IADR15(12), followed by the occupation
!              numbers in IADR15(13).
!              The calculation is performed for each of the NROOT
!              density matrices obtained in an average CASSCF calc.
!              Called from MAIN before OUTCTL
!
!          ****** IBM 3090 MOLCAS Release: 90 02 22 ******
!
      use rasscf_global, only: KSDFT, lRoots, NACPAR, NACPR2, iADR15,   &
     &                         iTri
      use SplitCas_Data, only: DoSPlitCas,lRootSplit
      use PrintLevel, only: DEBUG,USUAL
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,JOBIPH,NASH,NBAS,NFRO,NISH,NTOT,NTOT2
      Implicit None

      REAL*8 CMOO(*),SCR1(*),SCR2(*),SMAT(*),CMON(*),OCCN(*)

      Integer iPrLev, iDisk, jDisk, kRoot, I, IA, IB, IBAS, ID, IEND,   &
     &        II, IO, iOff, IST, iSTMO, ISYM, J, JA, jOff, NA1, NAO, NB,&
     &        NB2, NFI, ISTMO1

      IPRLEV=IPRLOC(7)
      iDisk=IADR15(12)
      jDisk=IADR15(3)

      if(.not.DoSplitCAS) then
        Do kRoot = 1,lRoots
          If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
            Write(LF,*)
            Write(LF,'(6X,A,I3)')                                       &
     &          'Natural orbitals and occupation numbers for root',kRoot
          End If
          Call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
          Call DBLOCK(SCR1)

          Call dCopy_(NTOT,[0.0d0],0,OCCN,1)
          Call dCopy_(NTOT2,CMOO,1,CMON,1)

          ID=0
          ISTMO1=0
          IB=0
          DO ISYM=1,NSYM
            NB=NBAS(ISYM)
            NFI=NFRO(ISYM)+NISH(ISYM)
            NAO=NASH(ISYM)
            NB2=NB**2
            NA1=ITRI(NAO+1)
            IO=IB+NFI
            ISTMO=ISTMO1+NB*NFI
!
!  set occupation number of frozen and inactive orbitals
!
            Call dCopy_(NFI,[2.0d0],0,OCCN(IB+1),1)
!
!  Diagonalize the density matrix and transform orbitals
!
            IF(NAO.GT.0) THEN
              Call dCopy_(NAO*NAO,[0.0d0],0,SCR2,1)
              Call dCopy_(NAO,[1.0d0],0,SCR2,NAO+1)
              CALL JACOB(SCR1(ID+1),SCR2,NAO,NAO)
              II=0
              DO I=1,NAO
                II=II+I
                OCCN(IO+I)=SCR1(II+ID)
              END DO
              IST=IO+1
              IEND=IO+NAO
              If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
                Write(LF,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')       &
     &                'sym',iSym,':',(OCCN(I),I=IST,IEND)
              End If
              CALL DGEMM_('N','N',                                      &
     &                    NB,NAO,NAO,                                   &
     &                    1.0d0,CMOO(ISTMO+1),NB,                       &
     &                    SCR2,NAO,                                     &
     &                    0.0d0,CMON(ISTMO+1),NB)
            END IF
!
            ID=ID+NA1
            ISTMO1=ISTMO1+NB2
            IB=IB+NB
          END DO
!
! ORTHOGONALIZE NEW MO'S
!
          CALL SUPSCH(SMAT,CMOO,CMON)
          CALL ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)
!
! Place NOs wrt decreasing order of OCCN
!***  GLM commented it off for this reordering is not consistent
!GLM throughout Molcas
!
!GLM          iOff=0
!GLM          jOff=0
!GLM          Do iSym=1,nSym
!GLM            NB=nBas(iSym)
!GLM            NAO=nAsh(iSym)
!GLM            IA=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
!GLM            JA=1+jOff+nFro(iSym)+nIsh(iSym)
!GLM            Call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
!GLM            iOff=iOff+NB**2
!GLM            jOff=jOff+NB
!GLM          End Do

       If(IPRLEV.ge.DEBUG) Then
           Write(LF,*)
           Write(LF,*) ' CMON in NATORB_RASSCF after ORDER_ARRAYS'
           Write(LF,*) ' ---------------------'
           Write(LF,*)
           ioff=0
           Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
              end do
              iOff = iOff + (iBas*iBas)
            end if
           End Do
       end if
!
! Write new molecular orbitals and occupation numbers to JOBIPH
!
!
          CALL DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
          CALL DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)
        End Do

      else   ! if DoSplitCAS (GLMJ)...
        If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
          Write(LF,*)
          Write(LF,'(6X,A,I3)')                                         &
     &     'Natural orbitals and occupation numbers for root',lRootSplit
        End If
        Call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
        Call DBLOCK(SCR1)

        Call dCopy_(NTOT,[0.0d0],0,OCCN,1)
        Call dCopy_(NTOT2,CMOO,1,CMON,1)

        ID=0
        ISTMO1=0
        IB=0
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          NFI=NFRO(ISYM)+NISH(ISYM)
          NAO=NASH(ISYM)
          NB2=NB**2
          NA1=ITRI(NAO+1)
          IO=IB+NFI
          ISTMO=ISTMO1+NB*NFI
!
!  set occupation number of frozen and inactive orbitals
!
          Call dCopy_(NFI,[2.0d0],0,OCCN(IB+1),1)
!
!  Diagonalize the density matrix and transform orbitals
!
          IF(NAO.GT.0) THEN
            Call dCopy_(NAO*NAO,[0.0d0],0,SCR2,1)
            Call dCopy_(NAO,[1.0d0],0,SCR2,NAO+1)
            CALL JACOB(SCR1(ID+1),SCR2,NAO,NAO)
            II=0
            DO I=1,NAO
              II=II+I
              OCCN(IO+I)=SCR1(II+ID)
            END DO
            IST=IO+1
            IEND=IO+NAO
            If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
               Write(LF,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')        &
     &               'sym',iSym,':',(OCCN(I),I=IST,IEND)
            End If
            CALL DGEMM_('N','N',                                        &
     &                  NB,NAO,NAO,                                     &
     &                  1.0d0,CMOO(ISTMO+1),NB,                         &
     &                  SCR2,NAO,                                       &
     &                  0.0d0,CMON(ISTMO+1),NB)
          END IF
!
         ID=ID+NA1
         ISTMO1=ISTMO1+NB2
         IB=IB+NB
       END DO
!
! ORTHOGONALIZE NEW MO'S
!
       CALL SUPSCH(SMAT,CMOO,CMON)
       CALL ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)
!
! Place NOs wrt decreasing order of OCCN
!
       iOff=0
       jOff=0
       Do iSym=1,nSym
          NB=nBas(iSym)
          NAO=nAsh(iSym)
          IA=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
          JA=1+jOff+nFro(iSym)+nIsh(iSym)
          Call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
          iOff=iOff+NB**2
          jOff=jOff+NB
       End Do
!
! Write new molecular orbitals and occupation numbers to JOBIPH
!
       CALL DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
       CALL DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)

      end if

      END Subroutine NatOrb_RASSCF
