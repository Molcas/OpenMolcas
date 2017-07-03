#ifndef _HAVE_EXTRA_

      Subroutine mnbrak(ax,bx,cx,fa,fb,fc,Error_for_t,rMP,xrMP,xxrMP,
     &                 xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,
     &                 nAtoms,nPert,Scratch_New,Scratch_Org,
     &                 iPrint_Errors)
      Implicit Real*8 (A-H,O-Z)
      Dimension EC(3,nij),AC(3,nij),C_o_C(3),R_ij(3)
      Dimension rMP(nij,0:nElem-1,0:nPert-1),xnrMP(nij,nElem)
      Dimension xrMP(nij,nElem),xxrMP(nij,nElem)
      Dimension Scratch_Org(nij*(2+lMax+1)),Scratch_New(nij*(2+lMax+1))
      External Error_for_t
      End Subroutine mnbrak

#endif
