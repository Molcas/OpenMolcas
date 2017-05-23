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
CPAM97      Subroutine Vibtrm(ndim,Vib1,Vib2,Umin,Umax,Teas,R,PotR,Title)
      Subroutine Vibtrm(ndim,Umin,Umax,Teas,R,PotR,Title)
C
C Calculation of transition matrix elements over vibrational
C wave functions of a given series of observables by a Simpson
C quadrature. The vibrational wave functions have been computed
C in vibrot and stored on units Vibwvs1 = 12 (potential 1) and
C Vibwvs2 = 13 (potential 2) for each rotational quantum number.
C The observables are read as input for a sequence of r-values
C and fitted to an analytical form in subroutine Pot.
C Integration is performed on a logarithmic scale (u=ln(r))
C between Umin and Umax with a grid size del corresponding to
C ndim integration steps. ndim has to be odd.
C Lifetimes are computed for the upper state (potential 2)
C no degeneracy factor is included in these values.
C Oscillator strengths are computed for all vibrational levels.
C
C     ********** MOLCAS Release 91 05 01 **********
C
      Implicit Real*8 (A-H,O-Z)
      Integer print
      Character*80 Title
      Character*1 hyph
#include "dimensions.fh"
      Parameter(nemax=nVib_Max,ndimmx=5000,
     &          lwork=2*nemax*ndimmx)
#include "intinp.fh"

#include "WrkSpc.fh"
#include "SysDef.fh"
      Dimension S(2*nemax**2+nemax),Obs(2*nemax**2+nemax)
C-POW Dimension nv1w(6),nv2w(6),Sw(6),
************************************************************************
      Dimension nv1w(2*nemax**2+nemax),nv2w(2*nemax**2+nemax),
     &          Sw(2*nemax**2+nemax),
     &          E(nemax,nemax),Taui(nemax,nemax),Tau(nemax),
     &          Osc(nemax,nemax)
      Dimension R(*),PotR(*)
      Character*80 TmpLine
*
      Call qEnter('VibTrm')
      ChkSum=0.0d0
C
      hyph='-'
CPAM97      write(6,1001) Title,Teas
CPAM971001  format(1h1,1x,'Matrix elements of observables'/1x,A80//
CPAM97     *       1x,'Asymtotic energy difference (au):',E14.6)
      Write(6,*)
      Call CollapseOutput(1,'Matrix elements of observable: '//Title)
      Write(6,'(a,e14.8)')' Asymptotic energy difference (au):',Teas
      print=2
c
c     check dimensions
c
      ne1=nvib1+1
      ne2=nvib21+1
      ndim1=ndim+1
      nwork=(ne1+ne2)*ndim1
      IERR=0
      if(nwork.gt.lwork) Then
        Write(6,*)'VIBTRM Error: NWork.gt.LWork'
        Write(6,'(1x,a,2i8)')'NWork,LWork:',NWork,LWork
        IERR=1
      end if
      if(ne1.gt.nemax) Then
        Write(6,*)'VIBTRM Error: NE1.gt.NEMax'
        Write(6,'(1x,a,2i8)')'NE1,NEMax:',NE1,NEMax
        IERR=1
      end if
      if(ne2.gt.nemax) Then
        Write(6,*)'VIBTRM Error: NE2.gt.NEMax'
        Write(6,'(1x,a,2i8)')'NE2,NEMax:',NE2,NEMax
        IERR=1
      end if
      if(NDim1.gt.NDimMx) Then
        Write(6,*)'VIBMAT Error: NDim1.gt.NDimMx'
        Write(6,'(1x,a,2i8)')'NDim1,NDimMx:',NDim1,NDimMx
        IERR=1
      end if
      If(IERR.ne.0) Then
        Call Abend
      End If

C Allocate memory
      Call GetMem('VecVib','Allo','Real',indexVib,lwork)
      Call GetMem('VecX','Allo','Real',indexX,ndim)

C Set up grid size
      del=(Umax-Umin)/(ndim-1)

C Loop over rotational quantum numbers
      JEndA=J1A
      JEndB=J1B
      If (iallrot.eq.1) Then
        JEndA=J2A
        JEndB=J2B
      End If
      Do 100 J1=J1A,JEndA
       Jad1=J1-J1A+1
      Do 100 J2=J1B,JEndB
       Jad2=J2-J1B+1
       Write(6,*)
       Write(TmpLine,1002) J1,J2
       Call CollapseOutput(1,TmpLine)
       Write(6,1003) (hyph,i=1,80)
1002   Format('Rotational quantum number for potential 1: ',I3,
     *        ', for potential 2: ',I3)
1003   Format(1x,80A1)

C Read vibrational functions for this pair of J-values.
       Do 35 ipot=1,2
        If(ipot.eq.1) then
         ne=ne1
         ist=1
         iadr1=iad12(Jad1)
         Do 20 nv=1,ne1
          Call DDafile(Vibwvs1,2,Work(indexVib+ist-1),ndim1,iadr1)
          E(1,nv)=Work(indexVib+ist+ndim-1)
          ist=ist+ndim1
20       Continue
        Else
         ne=ne2
         ist=ndim1*ne1+1
         iadr2=iad13(Jad2)
         Do 21 nv=1,ne2
          Call DDafile(Vibwvs2,2,Work(indexVib+ist-1),ndim1,iadr2)
          E(2,nv)=Work(indexVib+ist+ndim-1)+Teas
          ist=ist+ndim1
21       Continue
        Endif

C compute overlap matrix S
        ist1=-ndim1+ndim1*(ipot-1)*ne1
        nv12=(ipot-1)*(ne1**2+ne1)/2
        Do nv1=1,ne
         ist1=ist1+ndim1
         ist2=-ndim1+ndim1*(ipot-1)*ne1
         Do nv2=1,nv1
          ist2=ist2+ndim1
          nv12=nv12+1

C Set up scalar product and integrate
          Do i=1,ndim
           Work(indexX+i-1)=Work(indexVib+i+ist1-1)*
     &                      Work(indexVib+i+ist2-1)*R(i)**2
          End Do
          Call Simpsn(Work(indexX),del,ndim,S(nv12))
         End Do
        End Do

c check overlap matrix for non-orthogonality
       nv12=(ipot-1)*(ne1**2+ne1)/2
       Smax=0.0D0
       Do 32 nv1=1,ne
        Do 32 nv2=1,nv1
         nv12=nv12+1
         If(nv1.eq.nv2) go to 32
         If(abs(S(nv12)).gt.abs(Smax)) Smax=S(nv12)
32     Continue
       If(abs(Smax).gt.1.d-04) Write(6,1200) ipot,Smax
1200   Format(/1x,'*****Warning: non-orthogonality between vibrational',
     *                       1x,'wave functions for potential',I2
     *      /13x,'largest overlap matrix element is',d14.6)

C Print overlap matrix
       if ( Print.ge.1 ) then
          Write(6,*)
          Write(6,'(1x,A,I3)')
     &         'Overlap matrix for vibrational wave functions '//
     &         'for potential number',ipot
          iwr=0
          nv12=(ipot-1)*(ne1**2+ne1)/2
          Do nv1=1,ne
             Do nv2=1,nv1
                iwr=iwr+1
                nv12=nv12+1
                nv1w(iwr)=nv1-1
                nv2w(iwr)=nv2-1
                Sw(iwr)=S(nv12)
             End Do
          End Do
          Write(6,'(6(3X,2I3,F12.6))')
     &         (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
       End If
35    Continue
C
C     Compute overlap matrix between pot1 and pot2 functions
C
      nv11=0
      nv12s=(ne1**2+ne1+ne2**2+ne2)/2
      nv12=nv12s
      ist1=-ndim1
      Do 50 nv1=1,ne1
       ist1=ist1+ndim1
       nv11=nv11+nv1
       ist2=ndim1*(ne1-1)
       nv22=(ne1**2+ne1)/2
       Do 50 nv2=1,ne2
        ist2=ist2+ndim1
        nv12=nv12+1
        nv22=nv22+nv2
C       Set up scalar product and integrate
        Do 45 i=1,ndim
         Work(indexX+i-1)=Work(indexVib+i+ist1-1)*
     &                    Work(indexVib+i+ist2-1)*R(i)**2
45      Continue
        Call Simpsn(Work(indexX),del,ndim,S(nv12))
        S(nv12)=S(nv12)/sqrt(S(nv11)*S(nv22))
        ChkSum=ChkSum + S(nv12)
50     Continue
c
c      print overlap matrix
c
       Write(6,1310)
1310   Format(/1x,'Overlap matrix for pot-1 and pot-2 functions')
       nv12=nv12s-ne2
       Do 52 nv1=1,ne1
        nv12=nv12+ne2
        Write(6,'(6(3X,2I3,F12.6))')
     &       (nv1-1,i-1,S(i+nv12),i=1,ne2)
52     Continue
c
c      compute transition moment
c
       nv11=0
       nv12=0
       ist1=-ndim1
       Do nv1=1,ne1
        ist1=ist1+ndim1
        nv11=nv11+nv1
        ist2=ndim1*(ne1-1)
        nv22=(ne1**2+ne1)/2
        Do nv2=1,ne2
         ist2=ist2+ndim1
         nv12=nv12+1
         nv22=nv22+nv2
C        Set up scalar product and integrate
         do 55 i=1,ndim
          Work(indexX+i-1)=Work(indexVib+i+ist1-1)*
     &                     Work(indexVib+i+ist2-1)*PotR(i)*R(i)**2
55       Continue
         Call Simpsn(Work(indexX),del,ndim,Obsr)
         Obs(nv12)=Obsr/sqrt(S(nv11)*S(nv22))

C Taui is the contribution to the inverse lifetime from
C vibrational states nv1 (lower) and nv2 (upper)
         Taui(nv1,nv2)=21.419474D0*(E(2,nv2)-E(1,nv1))**3*Obs(nv12)**2
C computed oscillator strengths
         Osc(nv1,nv2)=2.D0*(E(2,nv2)-E(1,nv1))*Obs(nv12)**2/3.d0

cc         ChkSum=ChkSum + Taui(nv1,nv2)
        End Do
       End Do
c
c      write matrix elements
c
       Write(6,*)
       Write(6,'(1x,A)')
     &      'Transition moments over vibrational wave functions '//
     &      '(atomic units)'
       iwr=0
       Do nv1=1,ne1
          Do nv2=1,ne2
             iwr=iwr+1
             nv1w(iwr)=nv1-1
             nv2w(iwr)=nv2-1
             Sw(iwr)=Obs(iwr)
          End Do
       End Do
       Write(6,'(6(2X,I3,1X,I3,1X,F11.6))')
     &      (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
       Write(6,*)
       Write(6,'(1x,A)')
     &      'Energy differences for vibrational wave functions'//
     &      '(atomic units)'
       iwr=0
       Do nv1=1,ne1
          Do nv2=1,ne2
             iwr=iwr+1
             nv1w(iwr)=nv1-1
             nv2w(iwr)=nv2-1
             Sw(iwr)=E(2,nv2)-E(1,nv1)
          End Do
       End Do
       Write(6,'(6(2X,I3,1X,I3,1X,F11.6))')
     &      (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
        Write(6,'(1x,A)')
     &      'Energy differences for vibrational wave functions'//
     &      '(cm-1)'
       iwr=0
       Do nv1=1,ne1
          Do nv2=1,ne2
             iwr=iwr+1
             nv1w(iwr)=nv1-1
             nv2w(iwr)=nv2-1
             Sw(iwr)=(E(2,nv2)-E(1,nv1))*219474.51
          End Do
       End Do
      Write(6,'(6(2X,I3,1X,I3,1X,F11.1))')
     &      (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
C
C      Print oscillator strengths
C
       Write(6,'(1x,A)')
     &      'Oscillator strengths for vibrational wave functions'
       iwr=0
       Do nv1=1,ne1
          Do nv2=1,ne2
             iwr=iwr+1
             nv1w(iwr)=nv1-1
             nv2w(iwr)=nv2-1
             Sw(iwr)=Osc(nv1,nv2)
          End Do
       End Do
       Write(6,'(6(2X,I3,1X,I3,1X,E11.5))')
     &      (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
C
C      Print lifetime
C
       Write(6,*)
       Write(6,'(1x,A,/,A)')
     &      'Contributions to inverse lifetimes (ns-1)',
     &      'No degeneracy factor is included in these values.'
       iwr=0
       Do nv1=1,ne1
          Do nv2=1,ne2
             iwr=iwr+1
             nv1w(iwr)=nv1-1
             nv2w(iwr)=nv2-1
             Sw(iwr)=Taui(nv1,nv2)
          End Do
       End Do
       Write(6,'(6(2X,I3,1X,I3,1X,E11.5))')
     &      (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
C Compute lifetimes for upper state
       Write(6,1600)
1600   Format(/1x,'Lifetimes (in nano seconds)'/3x,'v',7x,'tau')
       Do 75 nv2=1,ne2
        Tau(nv2)=0.0D0
        Do 74 nv1=1,ne1
         Tau(nv2)=Tau(nv2)+Taui(nv1,nv2)
74      Continue
        Tau(nv2)=1.0D0/Tau(nv2)
cc        ChkSum=ChkSum + Tau(nv2)
        Write(6,1700) nv2-1,Tau(nv2)
1700    Format(1x,i3,f10.2)
75     Continue
       Call CollapseOutput(0,TmpLine)

C End of loop over rotational quantum number
100   Continue
      Call CollapseOutput(0,'Matrix elements of observable: '//Title)

C Deallocate memory
      Call GetMem('VecVib','Free','Real',indexVib,lwork)
      Call GetMem('VecX','Free','Real',indexX,ndim)

CPAM97CPAM97 New code begins:
CPAM97      Do J1=J1A,J2A
CPAM97      Do J2=J1B,J2B
CPAM97C Read vibrational function for pot-1, with rot q.n.=J1:
CPAM97       iadr1=iad12(J1-J1A+1)
CPAM97       Do nv=0,nvib1
CPAM97        Call Dafile(Vibwvs1,2,Vibbuf,RtoI*ndim,iadr1)
CPAM97        call dcopy_(ndim,Vibbuf,1,Vib1(1,nv),1)
CPAM97        ERoVib1(nv,J1)=Vibbuf(ndim+1)
CPAM97       Continue
CPAM97C Normalize vibrational function:
CPAM97       Do nv=0,nvib1
CPAM97         Do i=1,ndim
CPAM97           Vibbuf(i)=Vib1(i,nv)**2*R(i)**2
CPAM97         End Do
CPAM97         Call Simpsn(Vibbuf,del,ndim,ovlp)
CPAM97         xnrm=1.0d0/sqrt(ovlp)
CPAM97         Do i=1,ndim
CPAM97           Vib1(i,nv)=xnrm*Vib1(i,nv)
CPAM97         End Do
CPAM97       End Do
CPAM97C Check orthonormality:
CPAM97       ovlmax=0.0d0
CPAM97       Do nv1=1,nvib1
CPAM97         Do nv2=0,nv1-1
CPAM97           Do i=1,ndim
CPAM97             Vibbuf(i)=Vib1(i,nv1)*Vib1(i,nv2)*R(i)**2
CPAM97           End Do
CPAM97           Call Simpsn(Vibbuf,del,ndim,ovlp)
CPAM97           ovlmax=max(ovlp,ovlmax)
CPAM97         End Do
CPAM97       End Do
CPAM97       If(ovlmax.ge.1.0D-4) then
CPAM97         Write(6,*)' Warning from VIBTRM:'//
CPAM97     &     ' Vibrational wave functions for potential 1'
CPAM97         Write(6,*)' are not orthonormal. Max overlap:',ovlmax
CPAM97       End If
CPAM97C Read vibrational function for pot-2, with rot q.n.=J2:
CPAM97       iadr2=iad13(J2-J1B+1)
CPAM97       Do nv=0,nvib21
CPAM97        Call Dafile(Vibwvs2,2,Vibbuf,RtoI*ndim,iadr2)
CPAM97        call dcopy_(ndim,Vibbuf,1,Vib2(1,nv),1)
CPAM97        ERoVib2(nv,J2)=Vibbuf(ndim+1)
CPAM97       Continue
CPAM97C Normalize vibrational function:
CPAM97       Do nv=0,nvib21
CPAM97         Do i=1,ndim
CPAM97           Vibbuf(i)=Vib2(i,nv)**2*R(i)**2
CPAM97         End Do
CPAM97         Call Simpsn(Vibbuf,del,ndim,ovlp)
CPAM97         xnrm=1.0d0/sqrt(ovlp)
CPAM97         Do i=1,ndim
CPAM97           Vib2(i,nv)=xnrm*Vib2(i,nv)
CPAM97         End Do
CPAM97       End Do
CPAM97C Check orthonormality:
CPAM97       ovlmax=0.0d0
CPAM97       Do nv1=1,nvib21
CPAM97         Do nv2=0,nv1-1
CPAM97           Do i=1,ndim
CPAM97             Vibbuf(i)=Vib2(i,nv1)*Vib2(i,nv2)*R(i)**2
CPAM97           End Do
CPAM97           Call Simpsn(Vibbuf,del,ndim,ovlp)
CPAM97           ovlmax=max(ovlp,ovlmax)
CPAM97         End Do
CPAM97       End Do
CPAM97       If(ovlmax.ge.1.0D-4) then
CPAM97         Write(6,*)' Warning from VIBTRM:'//
CPAM97     &     ' Vibrational wave functions for potential 2'
CPAM97         Write(6,*)' are not orthonormal. Max overlap:',ovlmax
CPAM97       End If
CPAM97C Compute Franck-Condon overlaps and transition moment integrals:
CPAM97       Do nv1=0,nvib1
CPAM97         Do nv2=0,nvib21
CPAM97           If(J1.eq.J2) Then
CPAM97             Do i=1,ndim
CPAM97               Vibbuf(i)=Vib1(i,nv1)*Vib2(i,nv2)*R(i)**2
CPAM97             End Do
CPAM97             Call Simpsn(Vibbuf,del,ndim,FC(nv1,nv2,J1))
CPAM97           End If
CPAM97           If(abs(nv1-nv2).le.NSel) Then
CPAM97             Do i=1,ndim
CPAM97               X(i)=PotR(i)*Vibbuf(i)
CPAM97             End Do
CPAM97             Call Simpsn(X,del,ndim,Obs(nv1,nv2,J1,J2))
CPAM97           Else
CPAM97             TD(nv1,nv2,J1,J2)=0.0d0
CPAM97           End If
CPAM97           EDiff=ERoVib2(nv2,J2)-ERoVib1(nv1,J1)+TeDiff
CPAM97           TauPrt(nv1,nv2,J1,J2)=21.419474D0*EDiff**3*TD(nv1,nv2,J1,J2)**2
CPAM97         End Do
CPAM97       End Do
CPAM97      End Do
CPAM97      End Do
CPAM97CPAM97 New code ends
      Call Add_Info('VIBROT_VIBTRM',ChkSum,1,6)
      Call qExit('VibTrm')
      Return
      End
