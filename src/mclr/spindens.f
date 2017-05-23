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
      SubRoutine SpinDens
     &          (LS,RS,iL,iR,rP1,rp2,rp3,rp4,rp5,rDe1,rde2,
     &           itype)
*
      Implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "cstate_mclr.fh"
      Real*8 LS(*),RS(*),rP1(nna,nna,nna,nna),
     &       rP2(nna,nna,nna,nna),rP3(nna,nna,nna,nna),
     &       rP4(nna,nna,nna,nna),rP5(nna,nna,nna,nna),
     &       rDe1(nna,nna),rde2(nna,nna)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
*     Input:
*
*     LS : CI Coeff for left state
*     RS : CI Coeff for right state
*     iL : Symmetry of left state
*     iR : Symmetry of right state
*
*
*     Output:
*
*      rP1 : Two Electron -- Density
*      rP2 : Two Electron ++ Density
*      rP3 : Two Electron spin Density
*      rD1 : One Electron Spin adapted Density
*      rD1 : One Electron spin density Density
*
*
*     itype=1
*     Densities for:
*              0
*     <L|[{Q:S} ,H]|R>
*              0
*     itype 2
*     Densities for
*               0       0
*     <L||[{Q:S} ,[{Q:S} ,H]]|R>
*               0       0
*
C     Call GetMem('AAA','CHECK','REAL',i,i)

      Two=2.0d0
      rOne=1.0d0
      Zero=0.0d0
      n1=2*n1dens
      n2=2*n2dens+nnA**4
      Call GetMem('1Dens2','ALLO','Real',ipDe,n1)
      Call GetMem('2Dens2','ALLO','Real',ipP,n2)
      call dcopy_(n1,Zero,0,Work(ipDe),1)
      call dcopy_(n2,Zero,0,Work(ipP),1)
      nConfL=Max(ncsf(il),nint(xispsm(il,1)))
      nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))

      Call GetMem('CIL','ALLO','REAL',ipL,nConfL)
      Call GetMem('CIR','ALLO','REAL',ipR,nConfR)
      Call CSF2SD(LS,Work(ipL),iL)
      Call CSF2SD(RS,Work(ipR),iR)
      icsm=iR
      issm=iL
      Call Densi2(2,Work(ipDe),Work(ipP),
     &               Work(ipL),Work(ipR),0,0,1,n1dens,n2dens)
*
      If (itype.eq.1) Then
*
*             0
*        <0|[Q  ,H]|0>
*             0pq
*
         Do iA=1,nna
          Do jA=1,nnA
           rde1(ia,ja)=Work(ipDe+nna*(ia-1)+ja-1)-
     &                Work(ipDe+n1dens+nna*(ia-1)+ja-1)+
     &                Work(ipDe+nna*(ja-1)+ia-1)-
     &                Work(ipDe+n1dens+nna*(ja-1)+ia-1)
          End Do
         End Do
         Do iA=1,NnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             jilkab=(ja-1)*nna**3+(ia-1)*nna**2+(la-1)*nna+ka
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             jilkba=(la-1)*nna**3+(ka-1)*nna**2+(ja-1)*nna+ia
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             jilk=itri((ja-1)*nna+ia,(la-1)*nna+ka)
             rP1(ia,ja,ka,la)=Work(ipP+ijkl-1)-Work(ipP+ijkl-1+n2dens)+
     R        Work(ipP+n2dens*2+ijklab-1)-Work(ipP+ijklba-1+2*n2dens)+
     L        Work(ipP+jilk-1)-Work(ipP+jilk-1+n2dens)-
     L        Work(ipP+n2dens*2+jilkab-1)+Work(ipP+jilkba-1+2*n2dens)
            End Do
           End Do
          End Do
         End Do
*
*
      Else If (itype.eq.2) Then
* OK CONSTRUCT

*   --
*   ++
*   -+
*   spindensity
*   spin adadpted density

         Call DZAXPY(n1dens,-1.0d0,Work(ipDe),1,
     &           Work(ipDe+n1dens),1,rde1,1)
         Call DZAXPY(n1dens,1.0d0,Work(ipDe+n1dens),1,
     &           Work(ipDe),1,rde2,1)

         Do iA=1,nnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             ilkjab=(ia-1)*nna**3+(la-1)*nna**2+(ka-1)*nna+ja
             ilkjba=(ka-1)*nna**3+(ja-1)*nna**2+(ia-1)*nna+la
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             rP1(ia,ja, ka,la)=
     &                 Work(ipP+ijkl-1)+Work(ipP+n2dens+ijkl-1)-
     &                 Work(ipp+ijklab+2*n2dens-1)-
     &                 Work(ipp+ijklba+2*n2dens-1)
            End Do
           End Do
         End Do
         End Do
         Do iA=1,nnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             ilkjab=(ia-1)*nna**3+(la-1)*nna**2+(ka-1)*nna+ja
             ilkjba=(ka-1)*nna**3+(ja-1)*nna**2+(ia-1)*nna+la
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             rP2(ia,ja, ka,la)=
     &                 Work(ipP+ijkl-1)-Work(ipP+n2dens+ijkl-1)-
     &                 Work(ipp+ijklab+2*n2dens-1)+
     &                 Work(ipp+ijklba+2*n2dens-1)
            End Do
           End Do
         End Do
         End Do


         Do iA=1,NnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             rP3(ia,ja,ka,la)=Work(ipP+ijkl-1)+Work(ipP+ijkl-1+n2dens)+
     &        Work(ipP+n2dens*2+ijklab-1)+Work(ipP+ijklba-1+2*n2dens)
            End Do
           End Do
          End Do
         End Do
         Do iA=1,NnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             rP4(ia,ja,ka,la)=
     &        -Work(ipP+n2dens*2+ijklab-1)+Work(ipP+ijklba-1+2*n2dens)
            End Do
           End Do
          End Do
         End Do
         Do iA=1,NnA
          Do jA=1,nnA
           Do kA=1,nnA
            Do lA=1,nnA
             ijklab=(ia-1)*nna**3+(ja-1)*nna**2+(ka-1)*nna+la
             ijklba=(ka-1)*nna**3+(la-1)*nna**2+(ia-1)*nna+ja
             ijkl=itri((ia-1)*nna+ja,(ka-1)*nna+la)
             rP5(ia,ja,ka,la)=
     &        Work(ipP+n2dens*2+ijklab-1)+Work(ipP+ijklba-1+2*n2dens)
            End Do
           End Do
          End Do
         End Do
        End If
        Call GetMem('1Dens2','Free','Real',ipDe,n1)
        Call GetMem('2Dens2','Free','Real',ipP,n2)
        Call GetMem('CIL','Free','Real',ipL,nConfL)
        Call GetMem('CIR','Free','Real',ipR,nConfR)
      Return
      End
*
