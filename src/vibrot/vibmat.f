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
      Subroutine Vibmat(ne,ndim,Umin,Umax,Title,R,PotR,temp)
C
C     Calculation of matrix elements over vibrational wave functions
C     of a given series of observables by a simpson quadrature.
C     The vibrational wave functions have been computed in vibrot
C     and stored on unit Vibwvs=12 for each rotational quantum number.
C     The observables are read as input for a sequence of R-values
C     and fitted to an analytical form in function pot.
C     This routine is Called when the parameter iobs is not zero.
C     integration is performed on a logaritmic scale (u=ln(r))
C     between Umin and Umax with a grid size del corresponding to
C     ndim integration steps. ndim has to be odd.
C
C     ********** MOLCAS-2 Release 91 05 01 **********
C
      Implicit real*8 (a-h,o-z)
#include "constants.fh"
      Integer Print
      Real*8  Joule
      Logical WarnV,WarnR
      Character*80 Title
#include "dimensions.fh"
      Parameter(nemax=nVib_Max,lwork=nemax*500,
     &          ndimmx=5000,nem2=nemax*(nemax+1)/2)
#include "intinp.fh"
      Dimension Vib(lwork),S(nem2),X(ndimmx),
     *          Obs(nem2),R(*),PotR(*)
      Dimension nv1w(nem2),nv2w(nem2),Sw(nem2)
      Dimension SumO(nemax),SumW(nemax),Wmax(nemax),Wmin(nemax)
      Dimension Ener(nemax)
*
      Call qEnter('VibMat')
C
      Write(6,*)
      Call CollapseOutput(1,'Matrix elements of observable: '//Title)
      WarnV=.false.
      WarnR=.false.
      WarnVq=0.0d0
      WarnRq=0.0d0
      Boltz=CONST_BOLTZMANN_
      Joule=1.0d3*CONV_AU_TO_KJ_
      beta=1.0d0/(Boltz*Temp)
*     Write(6,'(a,e12.5)') 'Bolzmann ',Boltz
*     Write(6,'(a,e12.5)') 'au to J  ',Joule
*     Write(6,'(a,e12.5)') 'Temp     ',Temp
*     Write(6,'(a,e12.5)') 'beta     ',beta
      ndim1=ndim+1
      nwork=ne*ndim1
      IERR=0
      if(nwork.gt.lwork) Then
        Write(6,*)'VIBMAT Error: NWork.gt.LWork'
        Write(6,'(1x,a,2i8)')'NWork,LWork:',NWork,LWork
        IERR=1
      end if
      if(ne.gt.nemax) Then
        Write(6,*)'VIBMAT Error: NE.gt.NEMax'
        Write(6,'(1x,a,2i8)')'NE,NEMax:',NE,NEMax
        IERR=1
      end if
      if(NDim.gt.NDimMx) Then
        Write(6,*)'VIBMAT Error: NDim.gt.NDimMx'
        Write(6,'(1x,a,2i8)')'NDim,NDimMx:',NDim,NDimMx
        IERR=1
      end if
      If(IERR.ne.0) Then
        Call Abend
      End If

      Do nv=1,ne
         SumO(nv)=0.0d0
         SumW(nv)=0.0d0
         Wmax(nv)=0.0d0
         Wmin(nv)=0.0d0
      End Do

      del=(Umax-Umin)/(ndim-1)
      print=1
c
c     loop over rotational quantum numbers
c
      Do 100 J=J1A,J2A
C
C      Read vibrational functions for this J-value. store in vib.
C
       ist=1
       iad=iad12(J-J1A+1)
       Do nv=1,ne
        Call DDafile(Vibwvs,2,Vib(ist),ndim1,iad)
C       Write(6,6668) nv,(Vib(i+ist-1),i=1,ndim+1)
C6668   Format(/1x,'Vib-values',I3/(1x,6f12.6))
        ist=ist+ndim1
       Enddo
       Call DDafile(Vibwvs,2,Ener,ne,iad)
C
C      Compute overlap matrix S
C
       ist1=-ndim1
       nv12=0
       Do nv1=1,ne
        ist1=ist1+ndim1
        ist2=-ndim1
        Do nv2=1,nv1
         ist2=ist2+ndim1
         nv12=nv12+1
C
C        Set up scalar product and integrate
C
         Do i=1,ndim
          X(i)=Vib(i+ist1)*Vib(i+ist2)*R(i)**2
         Enddo
         Call Simpsn(X,del,ndim,S(nv12))
        Enddo
       Enddo
C
C      Check overlap matrix for non-orthogonality
C
       nv12=0
       Smax=0.0d0
       Do nv1=1,ne
        Do nv2=1,nv1
         nv12=nv12+1
         if(nv1.ne.nv2) Then
          if(abs(S(nv12)).gt.abs(Smax)) Smax=S(nv12)
         Endif
        Enddo
       Enddo
       if(abs(Smax).gt.1.d-04) write(6,1200) Smax
1200  Format(/1x,'***** Warning: Non-orthogonality between vibrational',
     *        1x,'wave functions.'
     *      /13x,'Largest overlap matrix element is',d14.6)
C
C      Print overlap matrix
C
       If ( Print.ge.2 ) then
       Write(6,'(1x,A)')
     & 'Overlap matrix for vibrational wave functions'
       iwr=0
       nv12=0
       Do nv1=1,ne
        Do nv2=1,nv1
         iwr=iwr+1
         nv12=nv12+1
         nv1w(iwr)=nv1
         nv2w(iwr)=nv2
         Sw(iwr)=S(nv12)
        End Do
       End Do
       Write(6,'(6(3X,2I3,F12.6))')
     & (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
       End If
c
c      set up integrand array for this Observable
c
       ist1=-ndim
       nv12=0
       nv11=0
       Do nv1=1,ne
        nv11=nv11+nv1
        nv22=0
        ist1=ist1+ndim
        ist2=-ndim
        Do nv2=1,nv1
         nv22=nv22+nv2
         ist2=ist2+ndim
         nv12=nv12+1
         Do i=1,ndim
          ri=R(i)
          X(i)=Vib(i+ist1)*Vib(i+ist2)*potR(i)*ri**2
         Enddo
         Call Simpsn(x,del,ndim,Obsr)
         Obs(nv12)=Obsr/sqrt(S(nv11)*S(nv22))
        Enddo
       Enddo
C
C      Write matrix elements
C
       Write(6,'(4h>>>>,1x,A,I3)')
     & 'matrix elements over vibrational wave functions '//
     & '(atomic units) for rotational quantum number',J
       iwr=0
       nv12=0
       Do nv1=1,ne
        Do nv2=1,nv1
         iwr=iwr+1
         nv12=nv12+1
         nv1w(iwr)=nv1
         nv2w(iwr)=nv2
         Sw(iwr)=Obs(nv12)
        Enddo
       Enddo
       Write(6,'(6(3X,2I3,F12.6))')
     & (nv1w(i),nv2w(i),Sw(i),i=1,iwr)
*     Dimension SumO(nemax),SumW(nemax),Wmax(nemax),Wmin(nemax)
*     Dimension Ener(nemax)
       Wfirst=1.0d0
       Wlast=0.0d0
       Do nv=1,ne
         ind=nv*(nv+1)/2
         O=Obs(ind)
         E=Ener(nv)*Joule
         W=exp(-beta*E)
         SumW(nv)=SumW(nv)+W
         SumO(nv)=SumO(nv)+W*O
         If(J.eq.J1a) Then
           Wmax(nv)=W
           Wmin(nv)=W
         Else
           Wmax(nv)=Max(W,Wmax(nv))
           Wmin(nv)=Min(W,Wmin(nv))
         End If
         If(nv.eq.1)  Wfirst=W
         If(nv.eq.ne) WLast=W
*        Write(6,'(a,3g15.6)') '... E,W,O ',E,W,O
       End Do
*      Write(6,'(a,f12.6)') 'Wlast/Wfirst',Wlast/Wfirst
       If(Wlast/Wfirst.gt.1.0e-3) Then
         WarnVq=Max(WarnVq,Wlast/Wfirst)
         WarnV=.true.
       End If
C
C     End of loop over rotational quantum number
C
100   Continue
      Do nv=1,ne
*        Write(6,'(a,i3)') 'Temp terms for vibration qn ',nv
*        Write(6,'(a,e15.6)') 'SumW(nv) ',SumW(nv)
*        Write(6,'(a,e15.6)') 'SumO(nv) ',SumO(nv)
*        Write(6,'(a,e15.6)') 'Wmax(nv) ',Wmax(nv)
*        Write(6,'(a,e15.6)') 'Wmin(nv) ',Wmin(nv)
*        Write(6,*)
         If(Wmin(nv)/Wmax(nv).gt.1.0d-3) Then
            WarnRq=Max(WarnRq,Wmin(nv)/Wmax(nv))
            WarnR=.true.
         End If
      End Do
      O=0.0d0
      W=0.0d0
      Do nv=1,ne
         O=O+SumO(nv)
         W=W+SumW(nv)
      End Do
      O=O/W
      Write(6,*)
      Write(6,'(a,e15.6,a,f9.3,a)') 'Temperature averaged observable:',
     &   O,'  at',Temp,'K'
      If(WarnV) Then
         Write(6,*)
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** Warning, temperature weighting not '//
     &                  'converged with respect to vibrational '//
     &                  'quantum numbers'
         Write(6,'(a)') '***'
         Write(6,'(a,e10.3,a)') '*** Quotient',WarnVq,' should be small'
         Write(6,'(a)') '***'
         Write(6,*)
      End If
      If(WarnR) Then
         Write(6,*)
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** Warning, temperature weighting not '//
     &                  'converged with respect to rotational '//
     &                  'quantum numbers'
         Write(6,'(a)') '***'
         Write(6,'(a,e10.3,a)') '*** Quotient',WarnRq,' should be small'
         Write(6,'(a)') '***'
         Write(6,*)
      End If
      Call CollapseOutput(0,'Matrix elements of observable: '//Title)
C
C     End of calculation for this Observable.
C
      Call qExit('VibMat')
      Return
      End
