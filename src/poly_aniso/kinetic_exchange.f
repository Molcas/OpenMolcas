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
      Subroutine Kinetic_Exchange(  N1,   N2,   M1,S1,  M2,S2,
     &                            eso1, eso2,
     &                            tpar, upar, lant, OPT, HKEX,
     &                            MR1,SR1,MR2,SR2)
c  compute KE, within various options :
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: lant,OPT
      Real(kind=8),intent(in)      :: tpar,upar
      !the Ln site
      Integer, intent(in)           :: N1
      Real(kind=8), intent(in)     :: eso1(N1)
      Complex(kind=8), intent(in)  :: M1(3,N1,N1)
      Complex(kind=8), intent(in)  :: S1(3,N1,N1)
      Complex(kind=8), intent(out) :: MR1(3,N1,N1)
      Complex(kind=8), intent(out) :: SR1(3,N1,N1)
      !the radical
      Integer, intent(in)           :: N2
      Real(kind=8), intent(in)     :: eso2(N2)
      Complex(kind=8), intent(in)  :: M2(3,N2,N2)
      Complex(kind=8), intent(in)  :: S2(3,N2,N2)
      Complex(kind=8), intent(out) :: MR2(3,N2,N2)
      Complex(kind=8), intent(out) :: SR2(3,N2,N2)
      ! exchange Hamiltonian
      Complex(kind=8) ::  HKEX(N1,N1,N2,N2)
      ! local variables:
      Integer          :: i1,i2,j1,j2,i,info,j,l,is1,is2,iprint
      Real(kind=8)    :: eloc1(N1)
      Real(kind=8)    :: eloc2(N2)
      Complex(kind=8) :: Z1(N1,N1)
      Complex(kind=8) :: Z2(N2,N2)
      Complex(kind=8) :: H1(N1,N1)
      Complex(kind=8) :: H2(N2,N2)
      Complex(kind=8) :: HCOV(N1,N1)
      Complex(kind=8) :: H1T(N1,N1)
      Complex(kind=8) :: HEXC( N1,N1,N2,N2)
      Complex(kind=8) :: ABIT( N1,N1,N2,N2)
      Complex(kind=8) :: MM1(3,N1,N1)
      Complex(kind=8) :: SM1(3,N1,N1)
      Complex(kind=8) :: ZZ1(N1,N1)
      Complex(kind=8) :: ZZ2(N2,N2)
      Complex(kind=8) :: ZCR(N1,N1)
      Complex(kind=8) :: TMP(N1,N1)
      Real(kind=8)    :: gtens(4,3),maxes(4,3,3),wcr(n1)
      Logical          :: DBG
      DBG=.false.
c determine the pseuDospin on each site (Z1 and Z2):
      Z1=(0.0_wp,0.0_wp)
      Z2=(0.0_wp,0.0_wp)
      iprint=1
      Call pseudospin(M1,N1,Z1,3,1,iprint)
      Call pseudospin(M2,N2,Z2,3,1,iprint)
      If (DBG) Then
        Call pa_prMat('KE_Exchange:: Pseudospin site 1',Z1,N1)
        Call pa_prMat('KE_Exchange:: Pseudospin site 2',Z2,N2)
      End If
c get the "traced-to-zero" local energy states on both sites:
      eloc1=0.0_wp
      eloc2=0.0_wp
      Call rtrace(N1,eso1,eloc1)
      Call rtrace(N2,eso2,eloc2)
      H1=(0.0_wp,0.0_wp)
      H2=(0.0_wp,0.0_wp)
      Do I1=1,N1
        Do I2=1,N1
           Do i=1,N1
           H1(I1,I2)=H1(I1,I2)+ELOC1(i)*conjg(Z1(i,I1))*Z1(i,I2)
           End Do
        End Do
      End Do
      Do I1=1,N2
        Do I2=1,N2
           Do i=1,N2
           H2(I1,I2)=H2(I1,I2)+ELOC2(i)*conjg(Z2(i,I1))*Z2(i,I2)
           End Do
        End Do
      End Do
      HCOV=(0.0_wp,0.0_wp)
      HEXC=(0.0_wp,0.0_wp)
      If ((OPT.eq.1).or.(OPT.eq.3).or.(OPT.gt.4)) Then ! full
      Call KE_Covalent(N1,   lant,tpar,upar,1,HCOV)
      Call KE_Exchange(N1,N2,lant,tpar,upar,1,HEXC)
      Else If ((OPT.eq.2).or.(OPT.eq.4)) Then ! 1/U model
      Call KE_Covalent(N1,   lant,tpar,upar,2,HCOV)
      Call KE_Exchange(N1,N2,lant,tpar,upar,2,HEXC)
      End If

      H1T=(0.0_wp,0.0_wp)
      H1T=H1+HCOV
c reWrite the HCOV in the initial ab initio basis:
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N1, N1,  N1, (1.0_wp,0.0_wp),
     &             Z1(1:N1,1:N1), N1,
     &           HCOV(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1 )
      HCOV=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N1, N1,  N1, (1.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1,
     &             Z1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &           HCOV(1:N1,1:N1), N1 )
c reWrite the H1 in the initial ab initio basis:
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &             Z1(1:N1,1:N1), N1,
     &             H1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1 )
      H1=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N1,  N1, N1, (1.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1,
     &             Z1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &             H1(1:N1,1:N1), N1 )
c reWrite the H2 in the initial ab initio basis:
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &             Z2(1:N2,1:N2), N2,
     &             H2(1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &            TMP(1:N2,1:N2), N2 )
      H2=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N2,  N2, N2, (1.0_wp,0.0_wp),
     &            TMP(1:N2,1:N2), N2,
     &             Z2(1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &             H2(1:N2,1:N2), N2 )
c reWrite the H1T in the initial ab initio basis:
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N1, N1,  N1, (1.0_wp,0.0_wp),
     &             Z1(1:N1,1:N1), N1,
     &            H1T(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1 )
      H1T=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N1, N1,  N1, (1.0_wp,0.0_wp),
     &            TMP(1:N1,1:N1), N1,
     &             Z1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &            H1T(1:N1,1:N1), N1 )
c reWrite the HEXC in the initial ab initio basis:
      Do is1=1,N2
        Do is2=1,N2
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &          Z1(1:N1,1:N1), N1,
     &   HEXC(1:N1,1:N1,is1,is2), N1, (0.0_wp,0.0_wp),
     &         TMP(1:N1,1:N1), N1 )
      HEXC(1:N1,1:N1,is1,is2)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N1,  N1, N1, (1.0_wp,0.0_wp),
     &         TMP(1:N1,1:N1), N1,
     &     Z1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &   HEXC(1:N1,1:N1,is1,is2), N1 )
        End Do
      End Do
      Do is1=1,N1
        Do is2=1,N1
      TMP(:,:)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &         Z2(1:N2,1:N2), N2,
     &   HEXC(is1,is2,1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &      TMP(1:N2,1:N2), N2 )
      HEXC(is1,is2,1:N2,1:N2)=(0.0_wp,0.0_wp)
      Call ZGEMM_('N','C',N2,  N2, N2, (1.0_wp,0.0_wp),
     &         TMP(1:N2,1:N2), N2,
     &     Z2(1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &     HEXC(is1,is2,1:N2,1:N2), N2 )
        End Do
      End Do
c
      HKEX=(0.0_wp,0.0_wp)
      ABIT=(0.0_wp,0.0_wp)
      Do i=1,N1
         Do i1=1,N2
      ABIT(i,i,i1,i1)=H1(i,i)+H2(i1,i1)
         End Do
      End Do
      Do i1=1,N1
        Do i2=1,N1
          Do j1=1,N2
            Do j2=1,N2
            HKEX(i1,i2,j1,j2)= HEXC(i1,i2,j1,j2) + ABIT(i1,i2,j1,j2)
            End Do
          End Do
        End Do
      End Do
c H1T = H1 + HCOV
      wcr=0.0_wp
      zcr=(0.0_wp,0.0_wp)
      info=0
      Call diag_c2(H1T,N1,info,wcr,zcr)
      Do i=1,N1
      Write(6,'(2(A,i2,A,F15.9))')
     & 'ESO1(',i,')=',ESO1(i),'  ESO1+COV(',i,')=',wcr(i)-wcr(1)
      End Do
      Do i=1,N1
       Do j=1,N1
      Write(6,'(A,i2,A,i2,A,2F20.14)') 'ZCR(',i,',',j,')=',ZCR(i,j)
       End Do
      End Do

c rotate to COV basis:
      If((opt.eq.3).or.(opt.eq.4)) Then
        Do is1=1,N2
          Do is2=1,N2
        TMP(:,:)=(0.0_wp,0.0_wp)
        Call ZGEMM_('C','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &           ZCR(1:N1,1:N1), N1,
     &    HKEX(1:N1,1:N1,is1,is2), N1, (0.0_wp,0.0_wp),
     &  TMP(1:N1,1:N1), N1 )
        HKEX(:,:,is1,is2)=(0.0_wp,0.0_wp)
        Call ZGEMM_('N','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &           TMP(1:N1,1:N1), N1,
     &     ZCR(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &   HKEX(1:N1,1:N1,is1,is2), N1 )
          End Do
        End Do
      End If
c compute the g tensors for initial and initial+covalence:
      Call atens(M1(1:3,1:2,1:2), 2, gtens(1,:), maxes(1,:,:), 1 )
      Call atens(M1(1:3,3:4,3:4), 2, gtens(2,:), maxes(2,:,:), 1 )
      MM1=(0.0_wp,0.0_wp)
      SM1=(0.0_wp,0.0_wp)
      Do L=1,3
        TMP(:,:)=(0.0_wp,0.0_wp)
        Call ZGEMM_('C','N',N1,N1,N1,(1.0_wp,0.0_wp), ZCR, N1,
     &           M1(L,:,:),N1,
     &             (0.0_wp,0.0_wp), TMP, N1 )
        Call ZGEMM_('N','N',N1,N1,N1,(1.0_wp,0.0_wp), TMP, N1, ZCR, N1,
     &             (0.0_wp,0.0_wp), MM1(L,:,:), N1 )
        TMP(:,:)=(0.0_wp,0.0_wp)
        Call ZGEMM_('C','N',N1,N1,N1,(1.0_wp,0.0_wp), ZCR, N1,
     &           S1(L,:,:),N1,
     &             (0.0_wp,0.0_wp), TMP, N1 )
        Call ZGEMM_('N','N',N1,N1,N1,(1.0_wp,0.0_wp), TMP, N1, ZCR, N1,
     &             (0.0_wp,0.0_wp), SM1(L,:,:), N1 )
      End Do
      Call atens(MM1(:,1:2,1:2), 2, gtens(3,:), maxes(3,:,:), 1 )
      Call atens(MM1(:,3:4,3:4), 2, gtens(4,:), maxes(4,:,:), 1 )
      Write(6,'(A)') 'Initial g tensors of the ground and first'//
     & ' excited KD'
      Do i=1,2
      Write(6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(i,1),
     & ' axis X: ',(maxes(i,j,1),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(i,2),
     & ' axis Y: ',(maxes(i,j,2),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(i,3),
     & ' axis Z: ',(maxes(i,j,3),j=1,3)
      Write(6,*)
      End Do
      Write(6,'(A)') 'Initial+Covalence g tensors of the ground and'//
     & ' first excited KD'
      Do i=3,4
      Write(6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(i,1),
     & ' axis X: ',(maxes(i,j,1),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(i,2),
     & ' axis Y: ',(maxes(i,j,2),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(i,3),
     & ' axis Z: ',(maxes(i,j,3),j=1,3)
      Write(6,*)
      End Do
c------
      MR1=(0.0_wp,0.0_wp)
      SR1=(0.0_wp,0.0_wp)
      MR2=(0.0_wp,0.0_wp)
      SR2=(0.0_wp,0.0_wp)
      If((opt.eq.3).or.(opt.eq.4)) Then
        ! rotate the magnetic moment to the coordinate
        ! system of main magnetic axes on Ln
        Call rotmom2(  SM1, N1, maxes(3,:,:), SR1 )
        Call rotmom2(  MM1, N1, maxes(3,:,:), MR1 )
        Call rotmom2(   S2, N2, maxes(3,:,:), SR2 )
        Call rotmom2(   M2, N2, maxes(3,:,:), MR2 )
        ! find the local pseuDospins on both sites:
        ZZ1=(0.0_wp,0.0_wp)
        ZZ2=(0.0_wp,0.0_wp)
        Call pseudospin(MR1,N1,ZZ1,3,1,iprint)
        Call pseudospin(MR2,N2,ZZ2,3,1,iprint)
        If (DBG) Then
          Call pa_prMat('KE_Exchange:: Pseudospin site 1',ZZ1,N1)
          Call pa_prMat('KE_Exchange:: Pseudospin site 2',ZZ2,N2)
        End If
!       rewrite the magnetic moments and spin moments in new local bases:
        Do L=1,3
          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N1,N1,N1,(1.0_wp,0.0_wp),
     &                ZZ1, N1,
     &                MR1(L,:,:), N1,  (0.0_wp,0.0_wp),
     &                TMP, N1 )
          MR1(L,:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('N','N',N1,N1,N1,(1.0_wp,0.0_wp),
     &               TMP, N1,
     &               ZZ1, N1, (0.0_wp,0.0_wp),
     &               MR1(L,:,:), N1 )

          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N1,N1,N1,(1.0_wp,0.0_wp),
     &               ZZ1, N1,
     &               SR1(L,:,:), N1, (0.0_wp,0.0_wp),
     &               TMP, N1 )
          SR1(L,:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('N','N',N1,N1,N1,(1.0_wp,0.0_wp),
     &               TMP, N1,
     &               ZZ1, N1, (0.0_wp,0.0_wp),
     &               SR1(L,:,:), N1 )
c
          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N2,N2,N2,(1.0_wp,0.0_wp),
     &               ZZ2, N2,
     &               MR2(L,:,:), N2, (0.0_wp,0.0_wp),
     &               TMP, N2 )
          MR2(L,:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('N','N',N2,N2,N2,(1.0_wp,0.0_wp),
     &               TMP, N2,
     &               ZZ2, N2, (0.0_wp,0.0_wp),
     &               MR2(L,:,:), N2 )
          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N2,N2,N2,(1.0_wp,0.0_wp),
     &               ZZ2, N2,
     &               SR2(L,:,:), N2, (0.0_wp,0.0_wp),
     &               TMP, N2 )
          SR2(L,:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('N','N',N2,N2,N2,(1.0_wp,0.0_wp),
     &               TMP, N2,
     &               ZZ2, N2, (0.0_wp,0.0_wp),
     &               SR2(L,:,:), N2 )
        End Do
c reWrite the exchnage matrix in the basis of local pseuDospins:
        Do is1=1,N2
          Do is2=1,N2
            TMP(:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('C','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &                  ZZ1(1:N1,1:N1), N1,
     &                  HKEX(1:N1,1:N1,is1,is2), N1, (0.0_wp,0.0_wp),
     &                  TMP(1:N1,1:N1), N1 )

            HKEX(1:N1,1:N1,is1,is2)=(0.0_wp,0.0_wp)
            Call ZGEMM_('N','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &                  TMP(1:N1,1:N1), N1,
     &                  ZZ1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &                  HKEX(1:N1,1:N1,is1,is2), N1 )
          End Do
        End Do
        Do is1=1,N1
          Do is2=1,N1
            TMP(:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('C','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &                  ZZ2(1:N2,1:N2), N2,
     &                  HKEX(is1,is2,1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &                  TMP(1:N2,1:N2), N2 )
            HKEX(is1,is2,1:N2,1:N2)=(0.0_wp,0.0_wp)
            Call ZGEMM_('N','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &                  TMP(1:N2,1:N2), N2,
     &                  ZZ2(1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &                  HKEX(is1,is2,1:N2,1:N2), N2 )
          End Do
        End Do
      Else !opt=1, opt=2, and opt>4
        MR1=M1
        SR1=S1
        MR2=M2
        SR2=S2
      End If
      Return
      End

