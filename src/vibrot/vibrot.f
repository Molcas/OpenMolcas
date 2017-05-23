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
      Subroutine Vibrot(ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,Redm,Req,
     &                  scale,temp)
C
C     Numerical solution of the vibrational Schroedinger equation
C     for a diatomic molecule.
C
C     PotR is the electronic energy given as input
C     Redm is the reduced mass
C     Jrot is the rotational quantum number
C     G is the radial wave function in a logaritmic scale U=ln(R)
C     Umin to Umax is the integration range
C     del is the step length
C     R=exp(U) is the radial coordinates
C     E is the energy
C     ndim is the number of integration steps
C     E0 is the starting value for an eigenvalue search. E0 is simply
C     put equal to the minimum value of the potential.
C     nvib is the number of eigenvalues wanted.
C     dE0 the original step length in the energy interpolation procedure
C
C     ********** MOLCAS-2 Release 90 05 01 **********
C
      Implicit Real*8 (A-H,O-Z)
C
#include "dimensions.fh"
#include "intinp.fh"
#include "observ.fh"
#include "warnings.fh"
      Dimension R(npoint+4),V(npoint),G(npoint),Svec(npoint),
     *          S(3*npoint),Ener(npoint),PotR(npoint+4),H(3*npoint),
     *          alfa(npoint),beta(npoint),W(npoint),B(npoint),
     *          Vec(npoint)
C
      Data zero/0.D 0/,one/1.D 0/,two/2.D 0/,three/3.D 0/,four/4.D 0/,
     *               five/5.D 0/,six/6.D 0/,twelve/12.D 0/
C
      D    = zero ! dummy initialize
      DMem = zero ! dummy initialize
C
      If(ngrid.gt.npoint) then
       Write(6,999) ngrid,npoint
999    Format(/1x,'Number of gridpoints',I4,' is larger than'/
     *        1x,' Dimension npoint=',I3,' in VIBROT.'/
     *        1x,' Program cannot continue.')
      Endif
C
C
      del=(Umax-Umin)/(ngrid-1)
      del2=del**2
      del12=del2/twelve
      del56=del2*five/six
      Thre=1.0D-08
      PotR0=PotR(ngrid+1)
      PotRn1=PotR(ngrid+2)
      R0=R(ngrid+1)
      Rn1=R(ngrid+2)
C
      Write(6,2700)
2700  Format(//1x,'Eigenstates'//1x,' Vib. q.n.   Rot. q.n.    Energy')
C
C     Start loop over rotational quantum numbers
C
      Est=E0
      Do 500 J=J1A,J2A
       iad12(J-J1A+1)=iadvib
       Xrot=-one/two+sqrt(one/four+J*(J+1)-lambda**2)
C
C      Construct effective potential V
C      and overlap vector Svec
C
       const=one/four+Xrot*(Xrot+1)
       redM2=two*redM
       Svec0=redM2*R0**2
       V0=const+Svec0*PotR0
       SvecN1=redM2*Rn1**2
       VN1=const+SvecN1*PotRn1
       Do 11 i=1,ngrid
        fact=redM2*R(i)**2
        Svec(i)=fact
        V(i)=const+fact*PotR(i)
11     Continue
C
C      Construct Hamilton and Overlap matrices
C
       k=3
       ngrid1=ngrid-1
       Do 20 i=2,ngrid1
        k=k+1
        H(k)=one-del12*V(i-1)
        S(k)=-del12*Svec(i-1)
        k=k+1
        H(k)=-two-DEL56*V(i)
        S(k)=-DEL56*Svec(i)
        k=k+1
        H(k)=one-del12*V(i+1)
        S(k)=-del12*Svec(i+1)
20     continue
       H(3)=one-del12*V(2)
       S(3)=-del12*Svec(2)
       nn1=3*ngrid-2
       H(nn1)=one-del12*V(ngrid1)
       S(nn1)=-del12*Svec(ngrid1)
C
C      H(1,1) and S(1,1) depends explicitly on energy and are
C      therefore constructed during the iteration process
C
C      set up and solve secular problem
C      search for eigenvalues starting from E0.
C
       fac=exp(-del*(Xrot+one/two))
       X2=redM/(two*Xrot+three)
       term2=fac*X2*R0**2
       term4=X2*R(1)**2
       H11a=-two-del56*V(1)
       Hnna=-two-del56*V(ngrid)
       H11c=one-del12*V0
       Hnnc=one-del12*VN1
       S11a=-del56*Svec(1)
       Snna=-del56*Svec(ngrid)
       S11c=-del12*Svec0
       Snnc=-del12*SvecN1
       delRN=Rn1-R(ngrid)
       nn=3*ngrid-1
C
       Do 100 ivib=1,nvib
        iter=0
        E=E0-dE0
        ist=0
        Emem=zero
        dE=dE0
C IFG Stop if above the dissociation limit
C     (assuming energy at infinity is set to zero)
101     If (E.gt.Zero) then
         Call SysQuitMsg(_RC_NOT_CONVERGED_,'VibRot',
     *        'Failed to find a state.','Try decreasing STEP.')
        End If
        E=E+dE
C       Compute h(1,1) and s(1,1)
        GAM0=(fac+(PotR0-E)*term2)/(one+(PotR0-E)*term4)
        H11=H11a+GAM0*H11c
        S11=S11a+GAM0*S11c
        H(2)=H11
        S(2)=S11
C       Compute H(N,N) and S(N,N)
        GAMn=exp(-sqrt(-redM2*E)*delRN-del/two)
        H(nn)=HnnA+GAMn*HnnC
        S(nn)=SnnA+GAMn*SnnC
C       Calculate secular determinant
        D1t=H11-E*S11
        D1=D1t/abs(D1t)
        D2=one
        k=3
        Do 110 i=2,ngrid
         Akk=H(k+2)-E*S(k+2)
         D=D1*Akk/abs(Akk)-D2*(H(k)-E*S(k))*(H(k+1)-E*S(k+1))/
     *   abs((H(k-1)-E*S(k-1))*Akk)
         D2=D1
         D1=D
         k=k+3
110     Continue
C
        If(ist.ne.0) go to 111
        ist=1
        Dmem=D
        go to 101
111     If(D/Dmem.lt.zero) go to 112
        Dmem=D
        go to 101
C       an eigenvalue has been bracketed
112     E1=E
        E=E1-DE*D/(D-Dmem)
        If(abs(E-Emem).lt.thre) go to 120
        Emem=E
        E=E-DE
        DE=DE/three
        E=E-DE
        ist=0
        go to 101
C       One eigenvalue found. Remove scaling
120     Ener(ivib)=(E-Est)/scale
        E0=E+dE0
        alfa(ngrid)=H(nn)-E*S(nn)
        beta(ngrid)=H(nn-1)-E*S(nn-1)
        k=nn+1
        Do 130 i=1,ngrid1
         k=k-3
         ii=ngrid-i
         If(ii.ne.1) beta(ii)=H(K-2)-E*S(K-2)
         Wii=(H(k)-E*S(k))/alfa(ii+1)
         alfa(ii)=H(k-1)-E*S(k-1)-Wii*beta(ii+1)
         W(ii)=Wii
130     Continue
C       step 2: unit right hand side
        If(abs(alfa(1)).eq.zero) alfa(1)=1.0D-10
        G(1)=one/alfa(1)
        Do 131 I=2,ngrid
         G(i)=(one-beta(i)*G(i-1))/alfa(i)
131     Continue
C       step 3: g on the right hand side
134     B(ngrid)=G(ngrid)
        Do 132 i=1,ngrid1
         ii=ngrid-i
         B(ii)=G(ii)-W(ii)*B(ii+1)
132     Continue
        Vec(1)=B(1)/alfa(1)
        Do 133 i=2,ngrid
         Vec(i)=(B(i)-beta(i)*Vec(i-1))/alfa(i)
133     Continue
C       normalize g and vec
        Gmax=zero
        Vmax=zero
        Do 140 i=1,ngrid
         AGi=abs(G(i))
         If(AGi.gt.Gmax) Gmax=AGi
         AVi=abs(Vec(i))
         If(AVi.gt.Vmax) Vmax=AVi
140     Continue
        diffm=zero
        Do 150 i=1,ngrid
         G(i)=G(i)/Gmax
         Vec(i)=Vec(i)/Vmax
         Adif=abs(Vec(i))-abs(G(i))
         If(Adif.gt.diffm) diffm=Adif
150     Continue
        If(diffm.lt.1.E-06.and.iter.ne.0) go to 160
        iter=iter+1
        Do 151 i=1,ngrid
         G(i)=Vec(i)
151     Continue
        go to 134
160     Continue
        Evib=(E-Est)/scale
CGG: Remove scaling
        Vec(ngrid+1)=E/scale
C        Vec(ngrid+1)=E
        Call DDafile(Vibwvs,1,Vec(1),ngrid+1,iadvib)
C
C       Print energies and notify on end point values
C
        V1=abs(Vec(1))
        Vn=abs(Vec(ngrid))
        THR=1.D-03
        If(V1.lt.THR.and.Vn.lt.THR) Write(6,3000) ivib-1,J,Evib
        If(V1.lt.THR.and.Vn.ge.THR) Write(6,3001) ivib-1,J,Evib,Vn
        If(V1.ge.THR.and.Vn.lt.THR) Write(6,3002) ivib-1,J,Evib,V1
        If(V1.ge.THR.and.Vn.ge.THR) Write(6,3003) ivib-1,J,Evib,V1,Vn
3000    Format(4X,I3,9X,I3,2X,F12.6)
3001    Format(4X,I3,9X,I3,2X,F12.6,5X,
     *         'Warning:function value at Rmax is',f12.6)
3002    Format(4X,I3,9X,I3,2X,F12.6,5X,
     *         'Warning:function value at Rmin is',f12.6)
3003    Format(4X,I3,9X,I3,2X,F12.6,5X,
     *         'Warning:function value at Rmin is',F12.6,' and at Rmax',
     *          F12.6)
C
100    Continue
C
       iadrsp(J-J1A+1)=iadvib
       Call DDafile(Vibwvs,1,Ener(1),nvib,iadvib)
       E0=Est
500   Continue
C
C     write address record to Vibwvs
C
      iadvib=0
      Call iDafile(Vibwvs,1,iad12,100,iadvib)
C
C     Compute spectroscopic parameters
C
      E0x=Est/scale
      If(ispc.ne.0) Call Spectc(Req,E0x,Atom1,Atom2,nvib)
C
C     Compute matrix elements of observables
C
      If(iobs.ne.0) then
       Do 600 i=1,iobs
       Call Vibmat(nvib,ngrid,Umin,Umax,Titobs(i),R,EoutO(1,i),temp)
600    Continue
      Endif
      Return
      End
