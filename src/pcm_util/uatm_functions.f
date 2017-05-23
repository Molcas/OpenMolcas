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
      Function AtNear(IAt,IAn,NBond,IBond,DH,DId,DAr,Chg)
      Implicit Real*8 (A-H,O-Z)
C     Logical AddH
      Parameter (MxBond=12)
      Dimension IAn(*),NBond(*),IBond(MxBond,*),Chg(*)
      AtNear= 0.0D+00
      NearAt= 0
      NH    = 0
      NNit  = 0
      NF    = 0
      NX    = 0
      NId   = 0
      NCSP3 = 0
      NCSP2 = 0
      NCSP  = 0
      DX    = 0.0D+00
      IAnI  = IAn(IAt)
      NbI   = NBond(IAt)
      Do 2000 j=1,NbI
        Jat = ibond(j,IAt)
        IAnJ= Ian(JAt)
        NbJ= NBond(JAt)
        ChgJ=Chg(JAt)
C Hydrogen
        If(IAnJ.eq.1) NH = NH + 1
C Carbon
        If(IAnJ.eq.6) then
          If(nbj.eq.4) NCSP3 = NCSP3 + 1
          If(nbj.eq.3) NCSP2 = NCSP2 + 1
          If(nbj.eq.2) NCSP  = NCSP  + 1
        EndIf
C Nitrogen
        If(IAnJ.eq.7) then
          IF(NBJ.lt.4) then
            If(IAnI.eq.IAnJ) NId = NId + 1
            NX=NX+1
            NNit=NNit+1
          else
            NX=NX+1
          endif
        endif
C Oxygen
        If(IAnJ.eq.8) then
          NOX = 1
          If(IAnI.eq.IAnJ.and.abs(ChgJ).lt.1.0D-01) NId = NId + 1
          If(CHGJ.gt.4.0D-01) then
            NOX = 0
            DX  = 5.0D-01
          EndIf
          If(CHGJ.lt.-4.0D-01) NOX = 3
          NX = NX + NOX
        EndIf
C Fluorine
        If(IAnJ.eq.9) then
          If(IAnI.eq.IAnJ) NId = NId + 1
          NF = NF + 1
        EndIf
C Positive Phosphorous
        If(IAnJ.eq.15.and.NBJ.eq.4) DX = DX - 1.5D0
C Positive Sulfur
        If(IAnJ.eq.16.and.NBJ.eq.3) DX = DX - 2.7D0
 2000 continue
      NearAt = NH
      NAr = 0
      If(IAn(IAt).eq.6.and.abs(Chg(IAt)).lt.4.0D-01)
     $  NearAt = NearAt + NH
      If(IAn(IAt).eq.6.and.NBond(IAt).eq.3.and.NCSP2.ge.1) then
        If((NH+NCSP3).ne.2)
     $    NAR = NAlpAr(IAt,IAn,NBond,IBond,Chg)
      EndIf
      If(IAn(IAt).eq.6.and.NBond(IAt).eq.4)
     $    NearAt =  NearAt - NX - NCSP2 - NNit + NF
     $                + NCAlph(IAt,NH,NCSP3,IAn,NBond,IBond,Chg)
      If(NH.gt.3) NearAt = NearAt - 1
      AtNear = DBLE(NearAt) - DX
      DH = DBLE(NH)
      DAr = 5.0D-01*DBLE(NAr)
      DId = 5.0D-01*DBLE(NId)
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Function Hybnew(OKUAH,OKCHG,IAt,IAn,NBond,IBond,IBType,PBO,Chg)
      Implicit Real*8 (A-H,O-Z)
      Logical PIAT,OKUAH,OKChg
      Parameter (MxBond=12)
      Dimension NBond(*),IBond(MxBond,*),IBtype(MxBond,*),IAn(*)
      Dimension PBO(MxBond,*)
C
      NumatI=IAn(IAt)
      HybNew = 3.0D+00
C     Chg    = 0.0D+00
      IQQ    = 0
      If(.not.OKUAH) then
        HybNew=0.0D+00
        return
      endif
C
C Hydrogen
C
      If(IAn(IAt).eq.1) then
        HybNew=0.0D+00
        If(NBond(IAT).eq.0.and..not.OKCHG) Chg = 1.0D+00
      EndIf
C
C Carbon
C
      If(IAN(IAt).eq.6) then
        NBTot=0
        PBtot=0.0D+00
        Do 2000 jj=1,NBond(IAt)
         NBtot=NBTot+IBType(jj,IAt)
         PBTot=PBTot+PBO(jj,IAt)
 2000 continue
        If(NBond(IAt).eq.3.and.NBTot.ge.4) HybNew= 2.0D+00
        If(NBond(IAt).eq.3.and.PBtot.gt.3.7d+0) Hybnew= 2.0d+00
        If(NBond(IAt).eq.2) HybNew = 1.0D+00
      EndIf
C
C Nitrogen, Phosphorous, Arsenicum, Antimonium
C
      If(IColAt(NumAtI).eq.5) then
        NBIA=NBond(IAt)
        goto(10,20,30,40)NBIA
  10    HybNew  = 1.0D+00
          ICC=IBond(1,IAt)
          ICN=IAn(IBond(1,IAt))
          If(ICN.eq.6.and.NBond(ICC).eq.1.and..not.OKCHG)
     $       Chg = -1.0D+00
          goto 50
  20    if(PIAT(IAt,IAn,NBond,IBond)) then
            HybNew = 2.0D+00
          else
            If(.not.OKCHG) Chg = -1.0D+00
          endif
          goto 50
  30    continue
          If(PIAT(IAt,IAn,NBond,IBond)) HybNew = 2.0D+00
          goto 50
  40    NBII=0
         Do 2010 ii=1,NBIA
         NumII=IAn(IBond(II,IAt))
         If(NumII.eq.6.or.NumII.eq.1) NBII=NBII+1
 2010 continue
        IF(NBII.gt.3.and..not.OKCHG) Chg = 1.0D+00
  50    continue
      EndIf
C
C Oxygen, Sulfur, Selenium, Tellurium
C
      If(IColAt(NumAtI).eq.6) then
C
C tricoordinated (positively charged only if bonded to C and/or H only)
C
        If(NBond(IAt).eq.3) then
         NBII=0
         do 2020 ii=1,3
           NumII=IAn(IBond(II,IAt))
           If(NumII.eq.1.or.NumII.eq.6) NBII=NBII+1
 2020 continue
         If(NBII.eq.3) then
           If(.not.OKCHG) Chg = 1.0D+00
           HybNew = 3.0D+00
         EndIf
        EndIf
C
C bicoordinated
C
        If(NBond(IAt).eq.2) then
          If(.not.OKCHG) Chg = 0.0D+00
          HybNew = 3.0D+00
C
C test for protonated aldehydes,ketones and similar for S
C
          Do 2030 jj=1,2
           JAt=IBond(JJ,IAt)
           NumJ=IAn(JAt)
           NBJ=NBond(JAt)
           If(NumJ.eq.6.and.NBJ.eq.3) then
             IQQ = 0
             do 2040 kk=1,3
                KAt = IBond(kk,JAt)
                NumK = IAn(KAt)
                NBK=NBond(KAt)
                If(NumK.eq.6.and.NBK.eq.4) IQQ=IQQ+1
 2040 continue
           endif
 2030 continue
          If(IQQ.ge.2) then
           If(.not.OKCHG) Chg = 1.0D+00
           HybNew = 2.0D+00
          EndIf
        EndIf
        If(NBond(IAt).eq.1) then
         HybNew = 2.0D+00
         Jat=IBond(1,IAt)
         icc=IAn(JAt)
         if(icc.eq.1) then
           HybNew = 3.0D+00
           If(.not.OKCHG) Chg = -1.0D+00
         EndIf
         if(icc.eq.6) then
           if(NBond(JAt).eq.4) then
             If(.not.OKCHG) CHg = -1.0D+00
             HybNew = 3.0D+00
           EndIf
           ISP2=0
           JSP2=0
           ICARB=0
           if(NBond(JAt).eq.3) JSP2=1
           do 2050 K=1,NBond(JAt)
             KAT=IBond(K,JAt)
             If(IAN(KAT).eq.6.and.NBond(KAT).eq.3) ISP2=ISP2+1
             IF(IAN(KAT).eq.8.and.NBond(KAt).eq.1) ICARB=ICARB+1
 2050 continue
           If(ISP2.gt.1) then
             If(.not.OKCHG) Chg = -1.0D+00
             HybNew = 3.0D+00
           EndIf
           If(JSp2.eq.1.and.ICARB.eq.2) then
             If(.not.OKCHG) Chg = -5.0D-01
             HybNew = 3.0D+00
           EndIf
         EndIf
        EndIf
      EndIf
C
C Halogens
C
      If(IColAt(NumAtI).eq.7.and.NBond(IAt).eq.0) then
        If(.not.OKCHG) Chg = -1.0D+00
      EndIf
C
C     Write(IOut,'(I3,2F3.1)') IAn(IAt),HybNew,Chg
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Integer Function NCAlph(IAt,NHI,NCSP3I,IAn,NBond,IBond,Chg)
      Implicit Real*8 (A-H,O-Z)
      Parameter (MxBond=12)
      Dimension NBond(*),IBond(MxBond,*),IAn(*),chg(*)
      NHetI=4-NHI-NCSP3I
      NCP1=0
      NCM1=0
C     write(IOut,'(''Atom'',I2,'' NHet='',I1)') IAt,NHetI
      Do 2000 jj=1,4
        jat=IBond(jj,iat)
        ianj=ian(jat)
        nbj=nbond(jat)
        nhj=0
        ncsp3j=0
        nhetj=0
        iplj=0
        if(ianj.eq.6.and.nbj.eq.4) then
          do 2010 kk=1,4
            kat=IBond(kk,JAt)
            iank=ian(KAt)
            nbk=nbond(kat)
            if(iank.eq.1) NHJ=NHJ+1
            if(iank.eq.6.and.nbk.eq.4) NCSP3J=NCSP3J+1
            if(chg(kat).gt.4.0d-01) iplj=1
 2010 continue
          NHETJ=4-NHJ-NCSP3J
          If(NHetI.ge.0.and.NHetJ.eq.0) NCP1=NCP1+1
          If(NHetI.eq.0.and.NHetJ.gt.0.and.iplj.eq.0) NCM1=NCM1+1
        endif
C       write(IOut,*) JAt,NHetJ,NCP1,NCM1
 2000 continue
      NCAlph = NCP1 - NCM1
      return
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Integer Function NAlpAr(IAt,IAn,NBond,IBond,Chg)
      Implicit Real*8 (A-H,O-Z)
      Parameter (MxBond=12)
      Dimension NBond(*),IBond(MxBond,*),IAn(*),Chg(*)
      NAlPar = -1
      NArI = 0
      Do 2000 jj=1,3
        NCSP2J=0
        NHetJ=0
        jat=IBond(jj,iat)
        ianj=ian(jat)
        nbj=nbond(jat)
        if(ianj.eq.7.and.nbj.ge.3) NCSP2J=2
        if(ianj.eq.6.and.nbj.eq.3) then
          do 2010 kk=1,3
            kat=IBond(kk,JAt)
            iank=ian(KAt)
            nbk=nbond(Kat)
            chgk=chg(KAt)
            if(chgk.lt.4.0D-01) then
             if(iank.eq.6.and.nbk.eq.3) NCSP2J=NCSP2J+1
             if(iank.eq.8.or.iank.eq.9) NHetJ=NHetJ+1
             if(iank.eq.17.or.iank.eq.35.or.iank.eq.53) NHetJ=NHetJ+1
             if(iank.eq.7) then
                if(nbk.le.2) NHetJ=NHetJ+1
                if(nbk.ge.3) NCSP2J=NCSP2J+1
             endif
            endif
 2010 continue
        endif
        if(NCSP2J.ge.2.and.NHetJ.eq.0) NarI = NArI + 1
 2000 continue
      If(NArI.ge.2) NAlPar = 1
      return
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Function IPBO(ToAng,IA,JA,RIJ,BondOr)
      Implicit Real*8(A-H,O-Z)
C
C     Return the order of the bond between atoms of atomic numbers IA
C     and JA, separated by RIJ, or 0 if there is no bond.
C
C     Generate connectivity based on bond distances alone.  The criteria
C     is whether the distances is no more than 30% longer than the sum of
C     covalent radii of involved atoms. For the moment all bond types are
C     determined using Pauling bond orders.
C     Note that RIJ is multiplied by ToAng, i.e. if it is already in Ang.
C     ToAng must be set = 1
C
      IPBO = 0
      R1IJ = RIJ*ToAng
      R0IJ = RCov97(IA,JA)
C     If(R1IJ.gt.R0IJ*1.3d0) return
      BondOr = Exp((R0IJ-R1IJ)/0.3d0)
      If(BondOr.lt.2.0D-01) return
      IBondO = Int(Sngl(BondOr+0.5d0))
      IBondO = Max(IBondO,1)
      IBondO = Min(IBondO,3)
      IPBO = IBondO
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Integer Function IColAt(NumAt)
      Dimension ICol(0:108)
      Save ICol
      Data ICol/0,1,2,
     $          1,2,3,4,5,6,7,8,
     $          1,2,3,4,5,6,7,8,
     $          1,2,9,9,9,9,9,9,9,9,9,9,3,4,5,6,7,8,
     $          1,2,9,9,9,9,9,9,9,9,9,9,3,4,5,6,7,8,
     $          1,2,30*9,
     $          22*9/
      IcolAt=ICol(NumAt)
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Logical Function PiAt(IAt,IAn,NBond,IBond)
      Parameter(MxBond=12)
      Dimension IAn(*),NBond(*),IBond(MxBond,*)
      PiAt=.False.
      ISPI=-1
      Do 2000 J=1,NBond(IAT)
        JAt=IBond(J,IAt)
        IAnJ=IAn(JAt)
        ITpJ=IColAt(IAnJ)
        NBnJ=NBond(JAT)
        NPiJ=0
        Do 2010 K=1,NBnJ
          KAt=IBond(K,JAt)
          If(IAn(KAt).eq.6.and.NBond(KAt).eq.3)NPiJ=NPiJ+1
 2010 continue
        If(IAnJ.eq.6.and.NBnJ.eq.3) then
          If(NPiJ.ge.2) then
            ISPI=ISPI+2
          Else
            ISPI=ISPI+1
          EndIf
        EndIf
        If(ITpJ.eq.5.and.NBnJ.eq.2) ISPI=ISPI+1
        If(ITpJ.eq.5.and.NPiJ.ge.2) ISPI=ISPI+1
 2000 continue
      If(ISPI.gt.0) PiAt=.True.
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Integer Function IRowAt(NumAt)
      Dimension IRow(0:108)
      Save IRow
      Data IROW /0,2*1,8*2,8*3,18*4,18*5,32*6,22*7/
      IRowAt=IRow(NumAt)
      Return
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Function AtSymb(I)
C                                    (I)
C     ATSYMB(I) is the Atomic Symbol corresponding to the Atomic number I
C     if I=0  Atsymb=Bq
C     if I=-1 Atsymb=X
      Character*2 AtSymb
#include "periodic_table.fh"
C
      If(I.gt.0) then
       AtSymb = PTab(i)
       Return
      EndIf
      If(I.eq.-1) AtSymb=' X'
      If(I.eq.0)  AtSymb='Bq'
      Return
      End
