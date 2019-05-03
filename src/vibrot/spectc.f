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
      Subroutine Spectc(Req,E0,Atom1,Atom2,nE)
C
C     Subroutine for calculation of spectroscopic constants for
C     diatomic molecules from the raw vibrational term values
C     obtained in Vibrot. A least square fit to these data is used
C
C     ********** MOLCAS-Release 91 05 01 **********
C
      Implicit Real*8 (A-H,O-Z)
#include "dimensions.fh"
      Parameter (ndim=nVib_Max*(nRot_Max+1),nemax=nVib_Max)
      Character*2 Atom1,Atom2
#include "intinp.fh"

#include "WrkSpc.fh"
#include "SysDef.fh"
C     Dimension Tvj(ndim),T(ndim),Gv(nemax),Gv2(nemax),F(100),Fc(100),
      Dimension Gv(nemax),Gv2(nemax),F(nemax),Fc(nemax),
     *          A(3,3),B(nemax),X(3),jind(nemax),D(nemax),G(3),dF(nemax)
*
      Call qEnter('Spectc')
C
      ChkSum=0.0d0
      zero=0.d0
      ntot=ne*(J2A-J1A+1)
      If(ntot.gt.ndim) then
        Write(6,*)'SPECTC Error: NTOT.gt.NDIM'
        Write(6,'(1x,a,2i6)')'NTOT,NDIM:',NTOT,NDIM
        Call Abend
      End If
      If(ne.gt.nemax) then
        Write(6,*)'SPECTC Error: NE.gt.NEMAX'
        Write(6,'(1x,a,2i6)')'NE,NEMAX:',NE,NEMAX
        Call Abend
      End If
c
      Call GetMem('VecTvj','Allo','Real',indexTvj,ndim)
      Call GetMem('VecT','Allo','Real',indexT,ndim)
c
      iend=0
      Do 10 J=J1A,J2A
       ist=iend+1
       iend=iend+ne
       Call DDafile(Vibwvs,2,Work(indexTvj+ist-1),ne,iadrsp(J-J1A+1))
10    Continue
C
C     Compute rotaional constants B and D for each vibrational band
C
      Const=219474.6D0
      If ( (J2A-J1A).lt.3 ) Then
        Write(6,*)'SPECTC Error: J2A-J1A.lt.3'
        Write(6,*)' SPECTC requires a range of rotational quanta'
        Write(6,*)' in order to fit spectroscopic constants.'
        Write(6,*)' The upper and lower limits of J must differ by'
        Write(6,*)'  at least 3.'
        Write(6,'(1x,a,2i6)')'J2A,J1A:',J2A,J1A
        Call Abend
      End If
      Do 100 nv=1,ne
      ij=ne
c     Term values for this band
      Do 90 j=J1A+1,J2A
       F(j)=Work(indexTvj+ij+nv-1)-Work(indexTvj+nv-1)
       F(j)=F(j)*const
       ij=ij+ne
90    Continue
      B1=0.D0
      B2=0.D0
      A11=0.D0
      A12=0.D0
      A22=0.D0
      Do 91 J=J1A+1,J2A
       Xj1=J*(J+1)-J1A*(J1A+1)
       Xj2=Xj1**2
       Xj3=Xj1**3
       Xj4=Xj2**2
       B1=B1+F(J)*Xj1
       B2=B2-F(J)*Xj2
       A11=A11+Xj2
       A12=A12+Xj3
       A22=A22+Xj4
91    Continue
c
      Xd=A12/A22
      Det=A11-Xd*A12
      B(nv)=(B1+B2*Xd)/Det
      D(nv)=(B2*A11/A22+B1*Xd)/Det
      Do 92 J=J1A+1,J2A
       Xj=J*(J+1)-J1A*(J1A+1)
       Fc(J)=B(nv)*Xj-D(nv)*Xj**2
       dF(J)=F(J)-Fc(J)
       nv1=nv-1
92    Continue
      write(6,990) nv1,B(nv),D(nv)
990   format(/5x,'Rotational constants for vibrational quantum number',
     *       i3/5x,'B=',e13.6,' cm-1     D=',e13.6,' cm-1'/
     *       5x,'Observed and computed term values (cm-1)')
      Do 93 J=J1A+1,J2A
       write(6,991) J,F(J),Fc(J),dF(J)
991    format(5x,i3,3e20.6)
93    Continue
100   Continue
      write(6,992)
992   format(/5x,'Rotational constants B(nv) and D(nv) in cm-1')
      Do 101 nv=1,ne
       write(6,991) nv,B(nv),D(nv)
101   Continue
C
C     Compute spectroscopic constants de (dele) and betae
C     by a least square fit
C
      dele=0.D0
      betae=0.D0
      if(ne.eq.1) de=d(1)
      if(ne.eq.1) go to 111
      a11=ne
      b1=0.D0
      b2=0.D0
      a12=0.D0
      a22=0.D0
      Do 110 nv=1,ne
       xnv=nv-0.5d0
       b1=b1+d(nv)
       b2=b2+d(nv)*xnv
       a12=a12+xnv
       a22=a22+xnv**2
110   Continue
      det=a11*a22-a12*a12
      dele=(b1*a22-b2*a12)/det
      betae=(b2*a11-b1*a12)/det
111   Continue
      write(6,994) dele,betae
994   format(/5x,'Spectroscopic constants De=',e13.6,' cm-1  Betae=',
     *       e13.6,' cm-1'/5x,'Observed and computed D values')
      Do 112 nv=1,ne
       nv1=nv-1
       dc=dele+betae*(nv-0.5d0)
       diff=d(nv)-dc
       write(6,991) nv1,d(nv),dc,diff
112   Continue
C
C     Spectroscopic constants Be,Alphe, and Gammae from rotational
C     constants B(nv)
C
      x(1)=0.D0
      x(2)=0.D0
      x(3)=0.D0
      if(ne.eq.1) x(1)=b(1)
      if(ne.eq.1) go to 126
      if(ne.eq.2) then
       x(1)=1.5D0*b(1)-0.5d0*b(2)
       x(2)=b(2)-b(1)
       go to 126
      Endif
      Do 115 i=1,3
       g(i)=0.D0
      Do 115 k=1,3
       a(i,k)=0.D0
115   Continue
      a(1,1)=ne
      Do 120 nv=1,ne
       xnv=nv-0.5d0
       xnv2=xnv**2
       xnv3=xnv**3
       xnv4=xnv2**2
       g(1)=g(1)+b(nv)
       g(2)=g(2)+b(nv)*xnv
       g(3)=g(3)+b(nv)*xnv2
       a(1,2)=a(1,2)+xnv
       a(1,3)=a(1,3)+xnv2
       a(2,3)=a(2,3)+xnv3
       a(3,3)=a(3,3)+xnv4
120   Continue
      a(2,1)=a(1,2)
      a(2,2)=a(1,3)
      a(3,1)=a(1,3)
      a(3,2)=a(2,3)
      call dminv(3,3,a)
      Do 125 i=1,3
       x(i)=0.D0
       Do 124 k=1,3
        x(i)=x(i)+a(i,k)*g(k)
124    Continue
125   Continue
126   Continue
      be=x(1)
      alphae=-x(2)
      gammae=x(3)
      write(6,995) be,alphae,gammae
995   format(/5x,'Spectroscopic constants Be,Alphae and Gammae'
     *       /5x,'Be=',e13.6,' cm-1    Alphae=',e13.6,' cm-1    Gammae='
     *       ,e13.6/5x,'Observed and computed B values')
      Do 127 nv=1,ne
       nv1=nv-1
       bc=be-alphae*(nv-0.5d0)+gammae*(nv-0.5d0)**2
       diff=b(nv)-bc
       write(6,991) nv1,b(nv),bc,diff
127   Continue
c
c     vibrational constants we,wexe and weye from band origins
c
      x(1)=0.D0
      x(2)=0.D0
      x(3)=0.D0
      if(ne.eq.1) x(1)=2*Work(indexTvj)
      if(ne.eq.1) go to 150
      if(ne.eq.2) then
       x(1)=3*Work(indexTvj)-Work(indexTvj+1)/3
       x(2)=2*Work(indexTvj+1)/3-2*Work(indexTvj)
       go to 150
      Endif
      Do 132 i=1,3
       g(i)=0.D0
      Do 132 k=1,3
       a(i,k)=0.D0
132   Continue
      Do 135 nv=1,ne
       xnv=nv-0.5d0
       xnv2=xnv**2
       xnv3=xnv**3
       g(1)=g(1)+Work(indexTvj+nv-1)*xnv
       g(2)=g(2)+Work(indexTvj+nv-1)*xnv2
       g(3)=g(3)+Work(indexTvj+nv-1)*xnv3
       a(1,1)=a(1,1)+xnv2
       a(1,2)=a(1,2)+xnv3
       a(1,3)=a(1,3)+xnv2**2
       a(2,3)=a(2,3)+xnv2*xnv3
       a(3,3)=a(3,3)+xnv3**2
135   Continue
      a(2,1)=a(1,2)
      a(2,2)=a(1,3)
      a(3,1)=a(1,3)
      a(3,2)=a(2,3)
      call dminv(3,3,a)
      Do 137 i=1,3
       x(i)=0.D0
       Do 136 k=1,3
        x(i)=x(i)+a(i,k)*g(k)
136    Continue
137   Continue
150   Continue
      we=x(1)*const
      wexe=x(2)*const
      weye=x(3)*const
      write(6,996) we,wexe,weye
996   format(/5x,'Vibrational constants we  =',e13.6,' cm-1'
     *       /5x,'                      wexe=',e13.6,' cm-1'
     *       /5x,'                      weye=',e13.6,' cm-1'
     *       /5x,'Observed and computed band origins')
      Do 155 nv=1,ne
       nv1=nv-1
       xnv=nv-0.5d0
       Tvjc=we*xnv+wexe*xnv**2+weye*xnv**3
       Tvjnv=const*Work(indexTvj+nv-1)
       diff=Tvjnv-Tvjc
       write(6,991) nv1,Tvjnv,Tvjc,diff
155   Continue
c
c     print output of spectroscopic constants
c
      write(6,1000) Atom1,Atom2
1000  format(//5x,'Spectroscopic constants for ',2a2)
      nv1=0
      nv2=ne-1
      write(6,1100) j1A,j2A,nv1,nv2
1100  format(//5x,'Range of J-values used in fit',2i3,
     *        /5x,'Range of v-values used in fit',2i3)
c
      Re=Req/1.88976d0
      dE=-E0*27.2099d0
      d0=dE-Work(indexTvj)*27.2099d0
c
      write(6,1200) Re,dE,d0,we,wexe,weye,
     *               Be,Alphae,Gammae,Dele,Betae
1200  format(///5x,'Re(a)',15x,f8.4,
     *         /5x,'De(ev)',14x,f8.4,
     *         /5x,'D0(ev)',14x,f8.4,
     *         /5x,'we(cm-1)',7x,e13.6,
     *         /5x,'wexe(cm-1)',5x,e13.6,
     *         /5x,'weye(cm-1)',5x,e13.6,
     *         /5x,'Be(cm-1)',7x,e13.6,
     *         /5x,'Alphae(cm-1)',3x,e13.6,
     *         /5x,'Gammae(cm-1)',3x,e13.6,
     *         /5x,'Dele(cm-1)',5x,e13.6,
     *         /5x,'Betae(cm-1)',4x,e13.6)
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7.>..+....8
      ChkSum=ChkSum + we + wexe
c
c     compute term values for check
c
      ind=0
      Do 300 J=J1A,J2A
       xj=J*(J+1)
       Do 290 nv=1,ne
        ind=ind+1
        xnv=nv-0.5d0
        Work(indexT+ind-1)=we*xnv+wexe*xnv**2+weye*xnv**3
     *        +xj*(be-alphae*xnv+gammae*xnv**2)
     *        -xj**2*(dele+betae*xnv)
290    Continue
300   Continue
c
c     compute max deviation
c
      tmax=zero
      Do 310 i=1,ntot
      Work(indexTvj+i-1)=Work(indexTvj+i-1)*const
      diff=abs(Work(indexT+i-1)-Work(indexTvj+i-1))
      if(diff.lt.tmax) go to 310
      tmax=diff
310   Continue
      ChkSum=ChkSum + tmax
      write(6,1300) tmax
1300  format(//5x,'Max deviation in term values is',e10.2,' cm(-1)')
c
c     print output of term values
c
      write(6,1400)
1400  format(//5x,'Term values(observed and computed) in cm(-1)')
      jst=j1A
315   jend=jst+4
      if(jend.gt.j2A) jend=j2A
      jrng=jend-jst+1
      Do 320 i=1,jrng
320   jind(i)=jst+i-1
      write(6,1500) (jind(i),i=1,jrng)
1500  format(/1x,'J-value',7x,i3,19x,i3,19x,i3,19x,i3,19x,i3)
      write(6,1510)
1510  format(/1x,'v-value')
      Do 340 nv=1,ne
      nv1=nv-1
      js=(jst-j1A)*ne+nv
      je=(jend-j1A)*ne+nv
      write(6,1600) nv1,(Work(indexTvj+i-1),Work(indexT+i-1),i=js,je,ne)
1600  format(4x,i2,2x,5(f8.2,2x,f8.2,4x))
340   Continue
      jst=jend+1
      if(jend.lt.j2A) go to 315
c
c     observed g-values (Tvj for j=0)
c
      Do 350 i=1,ne
      Gv(i)=Work(indexTvj+i-1)
      ChkSum=ChkSum + Gv(i)
350   Continue
c
c      compute delta g(v+1/2) values
c
      ne1=ne-1
      Do 420 i=1,ne1
420   Gv2(i)=Gv(i+1)-Gv(i)
c
c     print Gv and Gv2
c
      write(6,1700)
1700  format(//5x,'observed G-values in cm(-1)'/
     *       /5x,'v',6x,'G(v)',5x,'deltaG(v+1/2)')
      Do 440 i=1,ne1
      i1=i-1
      write(6,1800) i1,Gv(i)
1800  format(4x,i2,f11.2)
      write(6,1810) Gv2(i)
1810  format(17x,f13.2)
440   Continue
      i1=ne-1
      write(6,1800) i1,Gv(ne)
c
      Call GetMem('VecTvj','Free','Real',indexTvj,ndim)
      Call GetMem('VecT','Free','Real',indexT,ndim)
c
      Call Add_Info('VIBROT_SPECTC',[ChkSum],1,2)
*     Write(*,*) 'Spectc: ChkSum',ChkSum
      Call qExit('Spectc')
      Return
      End
