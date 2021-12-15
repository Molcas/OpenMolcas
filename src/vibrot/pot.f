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
      Subroutine POT(Rin,Ein,Rout,Eout,nout,ifit,Emin,Req,R0,R1,dR,
     *               npin,Title,iplot,Redm,scale,Nr)
      Implicit real*8 (A-H,O-Z)
      Parameter (lext=100,ipldim=1000)
      Character*80 Title
      Dimension Rout(*),Eout(*),Rext(lext),Eext(lext),Rin(*),Ein(*),
     *          iext(lext),Rplot(ipldim),Eplot(ipldim),Eextp(lext),
     *          Eplotp(ipldim)
#include "dimensions.fh"
#include "intinp.fh"
      Integer Print,Vibplt
      Character*15 type(3), FILENAME*8
      Data type/'Max point','Saddle point','Min point'/
      Data Vibplt/10/
C
*
      FILENAME='VIBPLT0 '
      If (Nr.ge.1 .and. Nr.lt.9) Then
         Write (FILENAME(7:7),'(I1)') Nr
      Else If (Nr.ge.10 .and. Nr.le. 99) Then
         Write (FILENAME(7:8),'(I2)') Nr
      End If
      Write (6,*) 'Generating plot file:',FILENAME
      Call Molcas_Open(Vibplt,FILENAME)
*
      Eeq=0.0D0  ! dummy intitialize
      Ree=0.0D0  ! dummy intitialize
*
      if(ifit.eq.1) goto 100
      if(ifit.eq.2) goto 200
      write(6,*)'POT Error: IFIT variable must be 1 or 2.'
      write(6,*)'           IFIT=',IFIT
      Call Quit_OnUserError()
100   Continue
      Ue=0.4D0
      print=1
      If(Ue.lt.0.1D0.or.Ue.gt.0.9D0) then
        write(6,*)'POT Error: Ue should be in 0.1..0.9'
        write(6,*)'           Ue  =',Ue
        Call Quit_OnUserError()
      end if
C
C     Find Re value
C
      Re=Rin(1)
      Eminx=Ein(1)
      Do 103 i=2,npin
       If(Ein(i).gt.Eminx) go to 103
       Re=Rin(i)
       Eminx=Ein(i)
103   Continue
C
C     Scale input potential if requested such that
C     the binding energy is 0.1 au (BOR 0601).
C
      scale=1.d0
      Escale=Ein(npin)
      if(iscale.ne.0) scale=abs(0.1d0/(Eminx-Escale))
      Do i=1,npin
       Ein(i)= Escale+scale*(Ein(i)-Escale)
      Enddo
      Redm=Redm/scale
      write(6,1002) scale
 1002 Format(/1x, 'Scaling parameter for potential:',f12.6)
C
      If(Re.lt.1.0D0.or.Re.gt.2.0D1) then
        write(6,*)'POT Error: Re should be in 1.0..20.0'
        write(6,*)'           Re  =',Re
        Call Quit_OnUserError()
      end if
      alpha=log(Ue)/Re
      Do i=1,npin
       Rin(i)=exp(alpha*Rin(i))
      Enddo
C
      Do  i=1,nout
        Rout(i)=exp(alpha*Rout(i))
      Enddo
      Rout(nout+1)=0.0D0
      next=0
      Call Sort_Pot(Rin,Ein,npin)
      next=lext
      Call Spline(Rin,Ein,npin,Rout,Eout,nout+1,Rext,Eext,iext,next,1)
      Do  i=1,npin
        Rin(i)=log(Rin(i))/alpha
      Enddo
      Call Sort_Pot(rin,ein,npin)
      Do i=1,nout
        Rout(i)=log(Rout(i))/alpha
      Enddo
      Do  i=1,next
        Rext(i)=log(Rext(i))/alpha
      Enddo
      If(print.ge.1) then
       Do i=1,next
        Eextp(i)=Escale+(Eext(i)-Escale)/scale
       Enddo
        if(next.ge.1) then
          Write(6,2002)
          Write(6,2003) (Type(iext(i)),Rext(i),
     *                   Eextp(i),i=1,next)
        End If
      End if
      If(iplot.ge.1) then
        nplot=1+int((R1-R0)/dR)
        If(nplot.gt.ipldim.or.nplot.le.0) then
          Write(6,*)'POT Error: Variable NPLOT should be in 1..IPLDIM'
          Write(6,*)'          IPLDIM=',IPLDIM
          Write(6,*)'          NPLOT =',NPLOT
          Call Quit_OnUserError()
        End If
        Do  i=1,nop
          Rin(i)=exp(alpha*Rin(i))
        Enddo
        Call Sort_Pot(Rin,Ein,nop)
        If(iplot.eq.1) then
          Do  i=1,nplot
            Rplot(i)=(i-1)*dR+R0
          Enddo
        Else
          Do  i=1,nplot
            Rplot(i)=1.0D1**((i-1)*dR+R0)
          Enddo
        End If
        Do  i=1,nplot
          Rplot(i)=exp(alpha*Rplot(i))
        Enddo
        next=lext
        Call Spline(Rin,Ein,npin,Rplot,Eplot,nplot,Rext,Eext,iext,
     *              next,1)
        Do  i=1,nplot
          Rplot(i)=log(Rplot(i))/alpha
        Enddo
        If(iplot.eq.3) then
          Do  i=1,nplot
            Rplot(i)=log10(Rplot(i))
          Enddo
        End If
        Do i=1,nplot
         Eplotp(i)=Escale+(Eplot(i)-Escale)/scale
        Enddo
        Write(Vibplt,3000) Title
        Write(Vibplt,3002) nplot
        Write(Vibplt,3001) (Rplot(i),Eplotp(i),i=1,nplot)
      End If
      if(next.le.0) then
        Write(6,*)'POT Error: Variable NEXT should be larger than 0'
        Write(6,*)'          NEXT =',NEXT
        Call Quit_OnUserError()
      End If
      Eeq=1.0d30
      imin=0
      Do  i=1,next
        if(Eeq.gt.Eext(i)) then
          imin=i
          Eeq=Eext(i)
        End If
      Enddo
      Emin=Eeq
      Ree=Rext(imin)
      if(iplot.gt.0) Ree=log(Ree)/alpha
      Req=Ree
      if(iext(imin).ne.3) then
        Write(6,*)'POT Error: IEXT(IMIN) should be = 3'
        Write(6,*)'    IEXT(IMIN) =',IEXT(IMIN)
        Call Quit_OnUserError()
      End If
      Einf=Eout(nout+1)
      Do  i=1,nout
        Eout(i)=Eout(i)-Einf
      Enddo
      Emin=Emin-Einf
      Einfp=Escale+(Einf-Escale)/scale
      Write(6,1950) Einfp
1950  Format(1x,'Extrapolated value at infinity',F13.6)
*
      Close(Vibplt)
*
      Return
C
C     Here for fitting of observable input
C
200   Continue
      Print=1
      Call Sort_Pot(Rin,Ein,npin)
      next=lext
      Call Spline(Rin,Ein,npin,Rout,Eout,nout,Rext,Eext,iext,next,1)
      If(iplot.ge.1) then
        nplot=1+int((R1-R0)/dR)
        If(nplot.gt.ipldim.or.nplot.le.0) then
          Write(6,*)'POT Error: Variable NPLOT should be in 1..IPLDIM'
          Write(6,*)'          IPLDIM=',IPLDIM
          Write(6,*)'          NPLOT =',NPLOT
          Call Quit_OnUserError()
        End If
        If(iplot.eq.1) then
          Do i=1,nplot
            Rplot(i)=(i-1)*dR+R0
          Enddo
        Else
          Do i=1,nplot
            Rplot(i)=1.0D1**((i-1)*dR+R0)
          Enddo
        End if
        next=lext
        Call Spline(Rin,Ein,npin,Rplot,Eplot,nplot,Rext,Eext,iext,
     *              next,1)
        If(iplot.eq.3) then
          Do i=1,nplot
            Rplot(i)=log10(Rplot(i))
          Enddo
        End if
        Write(Vibplt,3000) Title
        Write(Vibplt,3002) nplot
        Write(Vibplt,3001) (Rplot(i),Eplot(i),i=1,nplot)
      End if
      Req=Ree
      Emin=Eeq
*
      Close(Vibplt)
*
      Return
2002  Format(/1x,'extremum points'/24x,'R(au)',9x,'Value')
2003  Format(1x,A,2F14.6)
3000  Format(A80)
3001  Format(1x,f15.8,f20.8)
3002  Format(I4)
      END
