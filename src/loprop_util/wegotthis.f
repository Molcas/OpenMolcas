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
      Subroutine WeGotThis(nAt,nB,ipMP,ipC,nij,EC,iANr,ipT,ipTi,lMax
     &                    ,iTP,iPrint,Pot_Expo,Pot_Point,Pot_Fac,Diffed)
      Implicit real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension EC(3,nij),Pot_Expo(nij*2),Pot_Point(nij),Pot_Fac(nij*4)
      Dimension iANr(nAt)
      Dimension Expo(4),dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)

      Logical Diffed(nij*2)
      Logical Que,D1,D2

      Character*10 DistType(2),OneFile,Label

      Data DistType/'Monopole  ','Dipole    '/

*
*-- Print exponents and factors.
*
      Write(6,*)
      Write(6,'(A)')' *************************************************'
     &//'*********************************'
      Write(6,'(A)')' *                                                '
     &//'                                *'
      Write(6,'(A)')' *                  Result for the Diffuse Distrib'
     &//'ution                           *'
      Write(6,'(A)')' *                                                '
     &//'                                *'
      Write(6,'(A)')' *************************************************'
     &//'*********************************'
      Write(6,*)
      Write(6,401)'Centre','Coordinate','Multipole','Factor'
     &           ,'Exponent','Point-charge'
      Write(6,402)'|--------------------------------------------------'
     &          //'-----------------------------------------------|'
      kauntA=0
      Do iA=1,nAt
        Do jA=1,iA
          kauntA=kauntA+1
          Do iDC=1,2
            PP=-1.0d0
            If(iDC.eq.1)PP=Pot_Point(kauntA)
            If(Diffed(2*(kauntA-1)+iDC)) then
              If(iDC.eq.1) then
                Write(6,403)kauntA,(EC(k,kauntA),k=1,3),DistType(iDC)
     &                     ,Pot_Fac(4*(kauntA-1)+iDC)
     &                     ,2.0d0*Pot_Expo(2*(kauntA-1)+iDC),PP
              ElseIf(iDC.eq.2) then
                Write(6,404)kauntA,(EC(k,kauntA),k=1,3),DistType(iDC)
     &                     ,(Pot_Fac(4*(kauntA-1)+iDC+k),k=0,2)
     &                     ,2.0d0*Pot_Expo(2*(kauntA-1)+iDC)
              Endif
            Else
              If(iDC.eq.1) then
                Write(6,405)kauntA,(EC(k,kauntA),k=1,3),DistType(iDC)
     &                     ,Pot_Fac(4*(kauntA-1)+iDC),'Point(inf)',PP
              ElseIf(iDC.eq.2) then
                Write(6,406)kauntA,(EC(k,kauntA),k=1,3),DistType(iDC)
     &                     ,(Pot_Fac(4*(kauntA-1)+iDC+k),k=0,2)
     &                     ,'Point(inf)'
              Endif
            Endif
          Enddo
        Enddo
      Enddo
401   Format(' ',A,'     ',A,'             ',A,'             ',A
     &      ,'          ',A,'     ',A)
402   Format(A)
403   Format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,'          '
     &                            ,F7.3,'          ',F7.3,'      ',F7.3)
404   Format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,
     &              ' (',2(F7.3,','),F7.3,') ',F7.3)
405   Format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,'          '
     &                            ,F7.3,'            ',A,'      ',F7.3)
406   Format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,
     &              ' (',2(F7.3,','),F7.3,')   ',A)

*
*-- If the extra ONEINT-file exist, then do an error analysis.
*
      Que=.false.
      Call F_Inquire('ONEINTP',Que)
      If(Que) then
        Write(6,*)
        Write(6,'(A)')' Found Test-grid for error analysis.'
        Write(6,*)
        ErrAv1=0.0d0
        ErrAv2=0.0d0
        ErrDe1=0.0d0
        ErrDe2=0.0d0
        ErrVar1=0.0d0
        ErrVar2=0.0d0
        ErrMax1=0.0d0
        ErrMax2=0.0d0
        ErrCorr=0.0d0
        nImprove=0
        nShitty=0
        DeNom=0.0d0
        Write(OneFile,'(A)')'ONEINTP'
        Call Diff_Aux1(nEPP,ipEPCo,nB,OneFile)
        Call Get_D1ao(ip_D,nDens)
        Call GetMem('ElPot','Allo','Real',iElP,nDens+4)
        If(iPrint.ge.2) then
          Write(6,*)
          Write(6,'(A)')' Electric Potential'
          Write(6,'(A)')'  Reference   Approximate MP-expanded'
        Endif
*
*---- Loop over all points where the electric potential has
*     been sampled.
*
        Do iPP=1,nEPP
*
*---- First, get the true electric potential, the reference.
*
          Write(Label,'(A3,I5)')'EF0',iPP
          irc=-1
          iOpt=0
          iSmLbl=0
          iComp=1
          Call RdOne(irc,iOpt,Label,iComp,Work(iElP),iSmLbl)
          ElPot_REF=Work(iElP+nDens+3)
          ElPot_REF=ElPot_REF-Ddot_(nDens,Work(ip_D),1,Work(iElP),1)
*
*---- Second, get the approximate electric potential and also the
*     completely multipole expanded potential.
*
          ElPot_APP=0.0d0
          ElPot_MP=0.0d0
          kauntA=0
*          rMin=1d10
          Do iA=1,nAt
            Do jA=1,iA
              kauntA=kauntA+1
              x=Work(ipEPCo+(iPP-1)*3+0)-EC(1,kauntA)
              y=Work(ipEPCo+(iPP-1)*3+1)-EC(2,kauntA)
              z=Work(ipEPCo+(iPP-1)*3+2)-EC(3,kauntA)
              r=sqrt(x**2+y**2+z**2)
              rinv=1.0d0/r
              D1=Diffed(2*(kauntA-1)+1)
              D2=Diffed(2*(kauntA-1)+2)
              Expo(1)=Pot_Expo(2*(kauntA-1)+1)
              Expo(2)=Pot_Expo(2*(kauntA-1)+2)
              chP=Pot_Point(kauntA)
              kaunt=0
*              rmin=min(r,rmin)
              Do l=0,lMax
                kComp=(l+1)*(l+2)/2
                Do k=1,kComp
                  kaunt=kaunt+1
                  dMullig(kaunt)=Work(ipMP+nij*(kaunt-1)+kauntA-1)
                Enddo
              Enddo
              ElPot_APP=ElPot_APP+ElPot(r,rinv,x,y,z,dMullig,lMax,Expo
     &                                 ,chP,D1,D2)
              ElPot_MP=ElPot_MP+ElPot(r,rinv,x,y,z,dMullig,lMax,Expo
     &                                 ,chP,.false.,.false.)
            Enddo
          Enddo
*          Write(6,*)'Minimum Dist:',rMin

*
*---- Print if requested.
*
          If(iPrint.ge.2) then
            Write(6,441)ElPot_REF,ElPot_APP,ElPot_MP
     &                 ,(Work(ipEPCo+(iPP-1)*3+k),k=0,2)
          Endif
*
*---- Third, accumulate to error analysis.
*
*------ The difference
          Dif1=ElPot_APP-ElPot_REF
          Dif2=ElPot_MP-ElPot_REF
*------ Accumulate to average error
          ErrAv1=ErrAv1+Dif1
          ErrAv2=ErrAv2+Dif2
*------ Accumulate to variance of error
          ErrVar1=ErrVar1+Dif1**2
          ErrVar2=ErrVar2+Dif2**2
*------ Accumulate to deviation
          ErrDe1=ErrDe1+abs(Dif1)
          ErrDe2=ErrDe2+abs(Dif2)
*------ Maximal error
          If(ErrMax1.lt.abs(Dif1))ErrMax1=abs(Dif1)
          If(ErrMax2.lt.abs(Dif2))ErrMax2=abs(Dif2)
*------ Accumulate to covariance of errors
          ErrCorr=ErrCorr+Dif1*Dif2
*------ Better or worse
          If(abs(Dif1).le.abs(Dif2))nImprove=nImprove+1
          If(abs(Dif2).lt.abs(Dif1))nShitty=nShitty+1
*------ Accumulate to denominator in relative error
          DeNom=DeNom+abs(ElPot_REF)
        Enddo
        ErrRe1=ErrDe1/DeNom
        ErrRe2=ErrDe2/DeNom
        ErrAv1=ErrAv1/dble(nEPP)
        ErrAv2=ErrAv2/dble(nEPP)
        ErrDe1=ErrDe1/dble(nEPP)
        ErrDe2=ErrDe2/dble(nEPP)
        ErrVar1=ErrVar1/dble(nEPP)
        ErrVar2=ErrVar2/dble(nEPP)
        ErrCorr=ErrCorr/dble(nEPP)
        CorrCoeff=(ErrCorr-ErrAv1*ErrAv2)
     &           /sqrt((ErrVar1-ErrAv1**2)*(ErrVar2-ErrAv2**2))
        PImp=100.0*dble(nImprove)/dble(nEPP)
        PShi=100.0*dble(nShitty)/dble(nEPP)
*
*---- Four, print the analysis.
*
        Write(6,*)
        Write(6,'(A)')'   Error Analysis'
        Write(6,'(A)')' |----------------------------------------------'
     &//'------|'
        Write(6,'(A)')'                                Diffuse    '
     &//' MP-expanded'
        Write(6,442)'  Average absolute error:      ',ErrAv1,ErrAv2
        Write(6,442)'  Average absolute deviation:  ',ErrDe1,ErrDe2
        Write(6,442)'  Average relative error:      ',ErrRe1,ErrRe2
        Write(6,442)'  Maximal deviation:           ',ErrMax1,ErrMax2
        Write(6,*)
        Write(6,443)'  Error correlation:           ',CorrCoeff
        Write(6,444)'  Smaller error than MP:       ',PImp,'%'
        Write(6,444)'  Larger error than MP:        ',PShi,'%'
        Write(6,'(A)')' |----------------------------------------------'
     &//'------|'
        Write(6,*)
*
*---- Deallocate
*
        Call GetMem('ElPot','Free','Real',iElP,nDens+4)
        Call GetMem('Dens','Free','Real',ip_D,nDens)
        Call GetMem('PotPointCoord','Free','Real',ipEPCo,3*nEPP)
        irc=-1
        Call ClsOne(irc,0)
      Endif

441   Format(3(F12.8),'    In: ',3(F8.4))
442   Format(A,2(F12.8))
443   Format(A,F12.8)
444   Format(A,F12.8,A)

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(ipC)
         Call Unused_integer_array(iANr)
         Call Unused_integer(ipT)
         Call Unused_integer(ipTi)
         Call Unused_integer(iTP)
      End If
      End
