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
#include "stdalloc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "cstate_mclr.fh"
#include "real.fh"
      Real*8 LS(*),RS(*),rP1(nna,nna,nna,nna),
     &       rP2(nna,nna,nna,nna),rP3(nna,nna,nna,nna),
     &       rP4(nna,nna,nna,nna),rP5(nna,nna,nna,nna),
     &       rDe1(nna,nna),rde2(nna,nna)
      Real*8, Allocatable:: Dens(:,:), Pens(:), CIL(:), CIR(:)

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

      n1=2*n1dens
      n2=2*n2dens+nnA**4
      Call mma_allocate(Dens,n1dens,2,Label='Dens')
      Call mma_allocate(Pens,n2,Label='Pens')
      Dens(:,:)=Zero
      Pens(:)=Zero
      nConfL=Max(ncsf(il),NINT(xispsm(il,1)))
      nConfR=Max(ncsf(iR),NINT(xispsm(iR,1)))

      Call mma_allocate(CIL,nConfL,Label='CIL')
      Call mma_allocate(CIR,nConfR,Label='CIR')
      Call CSF2SD(LS,CIL,iL)
      Call CSF2SD(RS,CIR,iR)
      icsm=iR
      issm=iL
      Call Densi2(2,Dens,Pens,CIL,CIR,0,0,1,n1dens,n2dens)
*
      If (itype.eq.1) Then
*
*             0
*        <0|[Q  ,H]|0>
*             0pq
*
         Do iA=1,nna
          Do jA=1,nnA
           rde1(ia,ja)=Dens(nna*(ia-1)+ja,1)-
     &                 Dens(nna*(ia-1)+ja,2)+
     &                 Dens(nna*(ja-1)+ia,1)-
     &                 Dens(nna*(ja-1)+ia,2)
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
             rP1(ia,ja,ka,la)=Pens(ijkl)-Pens(ijkl+n2dens)+
     R        Pens(n2dens*2+ijklab)-Pens(ijklba+2*n2dens)+
     L        Pens(jilk)-Pens(jilk+n2dens)-
     L        Pens(n2dens*2+jilkab)+Pens(jilkba+2*n2dens)
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

         Call DZAXPY(n1dens,-One,Dens(:,1),1,Dens(:,2),1,rde1,1)
         Call DZAXPY(n1dens,1.0d0,Dens(:,2),1,Dens(:,1),1,rde2,1)

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
     &                 Pens(ijkl)+Pens(n2dens+ijkl)-
     &                 Pens(ijklab+2*n2dens)-
     &                 Pens(ijklba+2*n2dens)
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
     &                 Pens(ijkl)-Pens(n2dens+ijkl)-
     &                 Pens(ijklab+2*n2dens)+
     &                 Pens(ijklba+2*n2dens)
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
             rP3(ia,ja,ka,la)=Pens(ijkl)+Pens(ijkl+n2dens)+
     &        Pens(n2dens*2+ijklab)+Pens(ijklba+2*n2dens)
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
     &        -Pens(n2dens*2+ijklab)+Pens(ijklba+2*n2dens)
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
     &        Pens(n2dens*2+ijklab)+Pens(ijklba+2*n2dens)
            End Do
           End Do
          End Do
         End Do
        End If
        Call mma_deallocate(Dens)
        Call mma_deallocate(Pens)
        Call mma_deallocate(CIL)
        Call mma_deallocate(CIR)
      Return
      End
*
