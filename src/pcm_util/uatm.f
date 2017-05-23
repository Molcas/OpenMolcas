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
      Subroutine UATM(IOut,ICharg,NAt,NSfe,
     +                ToAng,Re,Alpha,C,IAn,NOrd,Chg,iPrint)
      Implicit real*8 (A-H,O-Z)
C
C New settings of radii for electrostatic cavity
C for HF/6-31(+)G* and ICOMP=4
C explicit values for C,N,O,F,S,Cl,Br,I, otherwise modified UFF radii
C
      Logical AlBond,OKUAH,OKCHG
C     Character AtSymb*2,AppNum*3,BCH*3,HH1*10,HH2*10,HY1*4,HY2*4,HY3*4
      Character AtSymb*2,         BCH*3,HH1*10,HH2*10,HY1*4,HY2*4,HY3*4
      Parameter (MxBond=12)

      Parameter (ISAX=1000)
      Dimension RE(*),C(3,NAt),IAn(*),NOrd(*),Chg(NAt)
      Dimension NBond(ISAX),IBond(MxBond,ISAX),IBtype(MxBond,ISAX)
      Dimension NTrBnd(ISAX),ITrBnd(MxBond,ISAX)
      Dimension ItrBtp(MxBond,ISAX)
      Dimension IHNum(ISAX),PBO(MxBond,ISAX)
      Dimension DHyb(ISAX),R0(0:7),Gamma(0:7),NOKUAH(20),Coeff(0:5)
      Dimension DX(ISAX),DAl(ISAX)
C     Save R0,Gamma,DRQM,DAlDon,DAlAcc,NOKUAH,Coeff
      Save R0,Gamma,DRQM,              NOKUAH,Coeff
      Data R0/0.0D+00,1.00D+00,1.50D+00,1.98D+00,2.08D+00,3*2.35D+00/
      Data Gamma/2*0.00D+00,9.00D-02,1.30D-01,4*1.50D-01/
C     Data DAlDon,DAlAcc,DRQM/2*1.0D-01,3.0D-01/
      Data               DRQM/          3.0D-01/
      Data NOKUAH/6,7,8,9,15,16,17,35,53,11*0/
      Data Coeff/1.0d+0,0.9d+0,0.6d+0,0.3d+0,0.1d+0,0.0d+0/
      BCH='sdt'
      HH1=' HHHHHHHHH'
      HH2='  23456789'
      HY1='*sss'
      HY2=' ppp'
      HY3='  23'
      Pi = 4.0d0*ATan(1.0d0)
C
C The number of explicitly parametrized atoms is NUAH and their atomic
C numbers are in NOKUAH
C
      NUAH=9
C
C find bonds
C
      OkChg = .False.
      AlBond=.False.
      call FndBnd(IOut,0,AlBond,ToAng,MxBond,NAt,
     +            IAn,C,NBond,IBond,IBtype,PBO,Re)
      NTotH=0
      Do 2010 IAt=1,NAt
       If(IAn(IAt).eq.1) NTotH=NTotH+1
       IHNum(IAt)=0
       NTrBnd(IAt)=0
       Do 2020 JJ=1,NBond(IAt)
        JHAT=IBond(JJ,IAT)
        If(IAn(JHAT).eq.1) then
         IHNum(IAt)=IHNum(IAt)+1
        Else
         NTrBnd(IAt)=NTrBnd(IAt)+1
         ITrBnd(NTrBnd(IAt),IAt)=JHAT
         ITrBtp(NTrBnd(IAt),IAt)=IBType(JJ,IAt)
        EndIf
 2020 continue
 2010 continue
      area=0.0d0
      do 2100 i1=1,nsfe
        iat=nord(i1)
        IEffbn=Min(5,NTrBnd(IAt))
        surf=4.0d0*Pi*re(I1)*re(I1)*coeff(IEffBn)
        area=area+surf
 2100 continue

      If (iPrint.gt.5) Then
         Write(iOut,*)
         Write(iOut,*)
         Write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
         Write(iOut,'(6X,A)') '================================'
         Write(iOut,*)
         write(IOut,'(6X,'' Nord Group'','//
     $   '''  Hybr  Charge Alpha Radius            Bonded to'')')
      End If
C
C Assign Charge and Hybridization to atoms
C
      QTot=0.0D+00
      Do 2030 IAt=1,NAt
       OKUAH=.False.
       do 2040 iverif=1,NUAH
        if(IAN(IAT).eq.NOKUAH(IVerif)) OKUAH=.true.
 2040 continue
       DHYB(IAt)=HybNew(OKUAH,OKCHG,IAt,IAn,NBond,IBond,IBtype,PBO,
     $  Chg(IAt))
       QTot=QTot+Chg(IAt)
 2030 continue
C
C Verify and possibly correct charges
C
      If (QTot.ge.1.0D-05) Then
         IQTot=Int(QTot+1.0D-01)
*     Else If(QTot.lt.1.0D-05)
      Else
         IQtot=Int(QTot-1.0D-01)
      End If
      If(IQTot.ne.ICharg) then
        FrQ = DBLE(ICharg)/DBLE(NAt-NTotH)
        Do 2050 IAt=1,NAt
C         write(IOut,'(6X,'' Atom'',I3,'' Old Q='',F5.2,
C    >       '' New Q='',F5.2)') IAt,Chg(IAt),FrQ
          If(IAn(IAt).ne.1) Chg(IAt) = FrQ
 2050 continue
      EndIf
C Determine Re, Alpha
      NSfe = 0
      Do 2060 IAt=1,NAt
        DDHYb = 0.0D+00
        DRQ   = 0.0D+00
        DDX   = 0.0D+00
        IRX  = IRowAt(IAn(IAt))
        RX   = R0(IRX)
        GX   = Gamma(IRX)
        NH   = 0
        NHI  = 0
        DH   = 0.0D+00
        DAr  = 0.0D+00
        DId  = 0.0D+00
        DX(IAt)= 0.0D+00
        OKUAH=.False.
        do 2070 iverif=1,NUAH
         if(IAN(IAT).eq.NOKUAH(IVerif)) then
           OKUAH=.true.
           DX(IAt)=
     $      AtNear(IAt,IAn,NBond,IBond,DH,DId,DAr,Chg)
         endif
 2070 continue
        DAl(IAt)= DAr - DId
C
C retain only H atoms not bonded to OKUAH atoms
C
        IUSE=1
        IF(IAn(IAt).eq.1) IUSE=0
C
        If(IUSE.eq.1) then
           NSfe = NSfe + 1
           I1 = NSfe
           NOrd(I1) = IAt
           OKUAH=.False.
             Re(I1)=Pauling(IAn(IAt))
             IAddH=Min(3,IHNum(IAt))
             Re(I1)=Re(I1)+DBLE(IAddH)*Gx
             do 2080 iverif=1,NUAH
              if(IAN(IAT).eq.NOKUAH(IVerif)) OKUAH=.true.
 2080 continue
           If(.not.okuah) goto 10
           If(Chg(IAt).lt.-1.0D-02) then
             DRQ = CHg(IAt)*DrQM
             IF(IAn(IAt).eq.7) DrQ = Chg(IAT)*2.0D-01
           EndIf
           If(Chg(IAt).gt.4.0D-01) then
             If(IAn(IAt).eq.8)  DRQ = -Chg(IAt)*2.6D-01
             If(IAn(IAt).eq.15) DRQ = -Chg(IAt)*4.5D-01
             If(IAn(IAt).eq.16) DRQ = -Chg(IAt)*5.5D-01
           EndIf
           NHI=Int(DHyb(IAt))
           NH=IHNum(IAt)
           If(NHI.lt.3) then
             DDHyb = (4.0D+00-DHyb(IAt))*GX
             If(IAn(IAt).ne.6) DDHyb=DDHyb/2.0D+00
           EndIf
           DDX = GX*(DX(IAT) + DAl(IAt))
           Re(I1) = RX + DDX + DDHyb + DRQ
C protect for too small C atoms
           IF(IAn(IAT).eq.6.and.Re(I1).lt.1.5D+00) Re(I1)=1.5D+00
  10       NH=IHNum(IAt)
             JE=Min(4,NTrBnd(IAt))
             If (iPrint.gt.5)
     &          Write(IOut,'(6X,1X,I3,2X,A2,2A1,3X,3A1,2X,'//
     $        'F5.2,3X,F4.2,2X,F5.3,3X,6(1X,A2,A3,''['',A1,'']''))')
     $        IAt,AtSymb(IAn(IAt)),HH1(NH+1:NH+1),HH2(NH+1:NH+1),
     $        HY1(NHI+1:NHI+1),HY2(NHI+1:NHI+1),HY3(NHI+1:NHI+1),
     $        Chg(IAt),Alpha,Re(I1),
     $        (AtSymb(IAn(ITrBnd(J,IAt))),'   ',
     $        BCH(ITrBtp(J,IAt):ITrBtp(J,IAt)),J=1,JE)
             If(NTrBnd(IAt).le.4) goto 20
             JE=Min(8,NTrBnd(IAt))
             If (iPrint.gt.5)
     &         Write(IOut,'(6X,40X,4(1X,A2,A3,''['',A1,'']''))')
     $         (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ',
     $         BCH(ITrBtp(JJ,IAt):ITrBtp(JJ,IAt)),JJ=5,JE)

             If(NTrBnd(IAt).le.8) goto 20
             JE=NTrBnd(IAt)
             If (iPrint.gt.5)
     &         Write(IOut,'(6X,40X,4(1X,A2,A3,''['',A1,'']''))')
     $         (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ',
     $         BCH(ITrBtp(JJ,IAt):ITrBtp(JJ,IAt)),JJ=9,JE)
  20         continue
        EndIf
 2060 continue
      If (iPrint.gt.5) Then
         write(IOut,'(6X,1X,78(''-''))')
         write(IOut,*)
      End If
      return
      end
