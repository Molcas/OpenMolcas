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
************************************************************************
      Subroutine LDKH(s,h,v,pvp,x,pxp,dkhorder,xorder,paramtype,
     &       dkhscfflg,N,isize,clight,vv,nn,dd,yy,ff,gg,xx,ii,jj,kk,mm,
     &       revt,eig,sinv,aa,rr,tt,pp,e,ew,snumber,tnumber,unumber,
     &       scrno1,scrno2,scr1,scr2,scr3,scr5,scr8,no_hamil,no_prop,
     &       nbasp,nbaso,indx2,Nsm,nblock,nAtom,scrsq,scrt,ssm,hsm,
     &       ewsm,eigsm)
************************************************************************
*                                                                      *
*     The local DKH procedure                                          *
*                                                                      *
*                                                                      *
*     Called from: dkh_driver                                          *
*                                                                      *
*     Calling    : sog                                                 *
*                  diag_dkh                                            *
*                                                                      *
*                                                                      *
*     Input      : h,s,v,pvp                                           *
*                                                                      *
*     Output     : h,x,pxp                                             *
*                                                                      *
************************************************************************
*      Implicit None
#include "dkhparameters.fh"
*
      Logical dkhscfflg,no_hamil,no_prop
      Integer  dkhorder,xorder,isize,snumber,tnumber,nbasp,nbaso,
     *        unumber,scrno1,scrno2,indx2(nAtom,4),Nsm,nblock,
     *        nAtom,N
      Real*8 clight,s(isize),h(isize),v(isize),pvp(isize),
     *                 x(isize),pxp(isize),ssm(Nsm,Nsm),
     *                 hsm(Nsm*(Nsm+1)/2)
      Real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      Real*8 revt(Nsm,Nsm),eig(Nsm,N),sinv(Nsm,Nsm),
     *                 aa(Nsm),rr(Nsm),tt(Nsm),pp(N),e(Nsm),
     *                 ew(N),ewsm(Nsm),eigsm(Nsm,Nsm)
      Real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(Nsm,N,2),scr8(isize,3),
     *                 scrsq(Nsm,Nsm,scrno1),
     *                 scrt(Nsm*(Nsm+1)/2,scrno2)
*
      Character*(3) paramtype
*
      Integer poss(maxsnumber),post(maxsnumber),posu(maxunumber)
*
      Integer i,j,k,l,m,counter,firstc,lastc,icen,jcen,itrian
*
      Real*8 det
      Integer adrmem,isfreeunit,dkh_48,adrnext
      Parameter (adrmem=4000)
      Integer dkhadr(adrmem)
#ifdef MOLPRO
      Character*32 filnam
#endif
      dkh_48=0

************************************************************************
*
      If (Out_Of_Core) then
        dkh_48=48
#ifdef _MOLCAS_
        dkh_48=isfreeunit(dkh_48)
        call daname_mf(dkh_48,'DKHMX')
#else
        filnam='DKHMX'
        call tmp_filename(filnam)
        call daname_nolib(dkh_48,filnam)
#endif
      End If
*
      If (dkhorder.gt.maxorder) then
        write(stdout,1010) dkhorder,maxorder,maxorder
1010    format (/,'Warning: dkhorder = ',I2,', but maxorder = ',I2,'.',
     *          /'--> dkhorder will be reduced to maxorder = ',I2,'.'/)
        dkhorder=maxorder
      End If
*
CMR150107      call initialize2 (N,isize,vv,nn,dd,yy,ff,gg,xx,ii,jj,kk,
CMR150107     *                  ll,mm,revt,tran,sinv,aa,rr,tt,pp,e,ew,snumber,
CMR150107     *                  tnumber,unumber,scrno1,scrno2,nbasp,
CMR150107     *                  scr1,scr2,scr3,scr5,scr8,dkhorder,xorder)
CMR260207 "nbaso" has not yet been included in initialize2
      Do i=1,maxsnumber
        poss(i)=0
        post(i)=0
      End Do
      Do i=1,maxunumber
        posu(i)=0
      End Do
*
************************************************************************
*
      Do i=1,isize
         scr8(i,1)=s(i)
      End Do
      Do i=1,isize
         s(i)=0.0d0
      End Do
      Do i=1,N
        Do j=1,N
          eig(i,j)=0.0d0
        End Do
      End Do
      icen=0
*                                                                      *
************************************************************************
*                 Create the small variable                            *
*                                                                      *
      Do i=1,nblock
        firstc=indx2(i,4)
        if (i.ne.nblock) Then
           lastc=indx2(i+1,4)-1
        Else
           lastc=nAtom
        End If
        counter=0
        itrian=1
        Do jcen=firstc,lastc
          counter=counter+indx2(jcen,2)-indx2(jcen,1)+1
          Do k=indx2(jcen,1),indx2(jcen,2)
             Do icen=firstc,jcen-1
                Do j=indx2(icen,1),indx2(icen,2)
                  ssm(itrian,1)=scr8(k*(k-1)/2+j,1)
                  hsm(itrian)=h(k*(k-1)/2+j)
                  itrian=itrian+1
                End Do
             End Do
             Do j=indx2(jcen,1),k
               ssm(itrian,1)=scr8(k*(k-1)/2+j,1)
               hsm(itrian)=h(k*(k-1)/2+j)
               itrian=itrian+1
             End Do
          End Do
        End Do
*
**    Verify that the overlap matrix is not singular
*
        call mat_sq_from_t (scrsq(1,1,1),counter,ssm)
        j=-1
        call mat_copy (scrsq(1,1,2),counter,counter,scrsq(1,1,1))
        call dcopiv (scrsq(1,1,2),scrsq(1,1,2),counter,1,counter,
     &               dkhzero,det,k,j,scrt(1,3))
        if (j.ne.0) then
          write (stdout,1020) dkhzero
1020      format('  D K H |****** '/,
     &  '        |****** WARNING - OVERLAP MATRIX SINGULAR '/,
     &  '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,
     &  '        |******'//)
          call errex_rel(' D K H | singular overlap matrix')
        endif
*
**      Run the calculation with the small variables
*
        Call sog(counter,ssm,sinv,scrt(1,1),scrt(1,2),ewsm)
        Call diag_dkh(hsm,counter,eigsm,ewsm,sinv,scrsq(1,1,2),1)
*
*      Reconstruct EIG and EW for the perturbation
*
        l=0
        Do icen=firstc,lastc
           Do k=indx2(icen,1),indx2(icen,2)
              ew(k)=ewsm(l+1)
              Do j=1,counter
                eig(j,k)=eigsm(l*counter+j,1)
                scr5(j,k,1)=sinv(l*counter+j,1)
                scr5(j,k,2)=scrsq(l*counter+j,1,1)
              End Do
              l=l+1
           End Do
        End Do
      End Do
*                                                                      *
************************************************************************
*                   Evaluation of the hamiltonian                      *
*
**      Reconstruct the small variables
*
      Do i=1,isize
         scr8(i,2)=v(i)
         scr8(i,3)=pvp(i)
         scr8(i,1)=x(i)
         s(i)=pxp(i)
         h(i)=h(i)+v(i)
      End Do
*
      Do i=1,nblock
        firstc=indx2(i,4)
        if (i.ne.nblock) Then
           lastc=indx2(i+1,4)-1
        Else
           lastc=nAtom
        End If
        counter=0
        l=0
        itrian=1
        Do icen=firstc,lastc
          counter=counter+indx2(icen,2)-indx2(icen,1)+1
        End Do
        Do icen=firstc,lastc
           Do k=indx2(icen,1),indx2(icen,2)
*
*             Faute ici aussi
*
              Do jcen=firstc,icen-1
                 Do j=indx2(jcen,1),indx2(jcen,2)
                   v(itrian)=scr8(k*(k-1)/2+j,2)
                   pvp(itrian)=scr8(k*(k-1)/2+j,3)
                   x(itrian)=scr8(k*(k-1)/2+j,1)
                   pxp(itrian)=s(k*(k-1)/2+j)
                   itrian=itrian+1
                 End Do
              End Do
              Do j=indx2(icen,1),k
                v(itrian)=scr8(k*(k-1)/2+j,2)
                pvp(itrian)=scr8(k*(k-1)/2+j,3)
                x(itrian)=scr8(k*(k-1)/2+j,1)
                pxp(itrian)=s(k*(k-1)/2+j)
                itrian=itrian+1
              End Do
              ewsm(l+1)=ew(k)
              Do j=1,counter
                eigsm(l*counter+j,1)=eig(j,k)
                sinv(l*counter+j,1)=scr5(j,k,1)
                ssm(l*counter+j,1)=scr5(j,k,2)
              End Do
              l=l+1
           End Do
        End Do
        isize=counter*(counter+1)/2
        if (.not.Out_Of_Core) nbaso=counter
*
**    Evaluation of the hamiltonian
*
        Call calc_prefactors (counter,isize,clight,aa,rr,tt,pp,e,ewsm)
        Call calc_revt (counter,revt,eigsm,sinv,ssm,scrsq(1,1,2))
        Call calc_E0 (counter,isize,hsm,revt,ewsm)
        Call calc_E1 (counter,isize,xorder,dkhscfflg,hsm,v,pvp,x,pxp,
     &                vv,dd,xx,jj,revt,eigsm,sinv,aa,rr,tt,scrno1,
     &                scrno2,scrsq,scrt,no_prop,nbasp,nbaso,dkhadr,
     &                adrmem,dkh_48,adrnext)
*       nn,yy,ii and kk are respectively vv,dd,xx and jj
       Call setup_matrices (counter,isize,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,
     &                     kk,ll,mm,aa,rr,e,scrno2,scrt,no_prop,
     &                     nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
       Call calc_Uxxx (counter,dkhorder,xorder,dkhscfflg,posu,post,poss,
     &                vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     &                tnumber,unumber,scrno1,scr1,scr2,scr3,scrsq,nbasp,
     &                nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
       Call calc_Txxx (counter,dkhorder,xorder,dkhscfflg,posu,post,poss,
     &                vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     &                tnumber,unumber,scrno1,scr1,scr2,scr3,scrsq,nbasp,
     &                nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
       Call calc_Sxxx (counter,dkhorder,xorder,dkhscfflg,posu,post,poss,
     &                vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     &                tnumber,unumber,scrno1,scr1,scr2,scr3,scrsq,nbasp,
     &                nbaso,dkhadr,adrmem,dkh_48,adrnext)

        If (.not.no_hamil)
     & Call calc_operators(counter,isize,dkhorder,xorder,dkhscfflg,posu,
     &                      post,poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,
     &                      kk,ll,mm,e,snumber,tnumber,unumber,
     &                      scrno1,scrno2,scr1,scr2,scr3,scrsq,scrt,
     &                      revt,hsm,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     &                      adrnext)
        If (.not.dkhscfflg.and..not.no_prop)
     &  Call calc_xoperators(counter,isize,dkhorder,xorder,
     &                       dkhscfflg,posu,post,poss,vv,nn,dd,yy,
     &                       ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     &                       tnumber,unumber,scrno1,scrno2,scr1,
     &                       scr2,scr3,scrsq,scrt,revt,x,nbasp,nbaso,
     &                       dkhadr,adrmem,dkh_48,adrnext)
*
**     Reconstruct the full variables
*
        itrian=1
        m=1
        l=1
        Do icen=firstc,lastc
           Do k=indx2(icen,1),indx2(icen,2)
              Do jcen=firstc,icen-1
                Do j=indx2(jcen,1),indx2(jcen,2)
                  h(k*(k-1)/2+j)=hsm(itrian)
                  itrian=itrian+1
                End Do
              End Do
              Do j=indx2(icen,1),k
                h(k*(k-1)/2+j)=hsm(itrian)
                itrian=itrian+1
              End Do
           End Do
        End Do
      End Do
*
      isize=N*(N+1)/2
      if (Out_Of_Core) then
#ifdef _MOLCAS_
        call DaEras(dkh_48)
#else
        call daclos(dkh_48)
#endif
      endif
*
c Avoid unused argument warnings
      If (.False.) Call Unused_character(paramtype)
      End
************************************************************************
************************************************************************
************************************************************************
      subroutine calc_indx(indx2,indx,Coord,N,nAtom,Nsm,nblock)
************************************************************************
*                                                                      *
*     Create the auxiliary index                                       *
*                                                                      *
*     Called from : dkrelint                                           *
*                                                                      *
*     Calling     : none                                               *
*                                                                      *
*                                                                      *
*     indx2( ,3) is the number of the block. It is equal to the number *
*     of the atom, or to the number of the center define by the user   *
*     when the atom is inside the radius                               *
*     indx2( ,1) and indx2( ,2) are respectively the number of the     *
*     first and the last primitive basis                               *
*     indx( ,4) is a kind of new array of the size nblock, which       *
*     contains the number of atoms in each block                       *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
      use DKH_Info
      Integer indx2(nAtom,4),indx(N),N,Nsm,nblock,itmp1,itmp2
      Real*8 Coord(3*nAtom),distce
*
      call get_iarray('Ctr Index Prim',indx,N)
      indx2(1,1)=1
      indx2(1,3)=1
      Do j=1,N
         indx2(indx(j),2)=j
      End Do
      Do j=2,nAtom
         indx2(j,3)=j
         indx2(j,1)=indx2(j-1,2)+1
      End Do
*
**    Make one block with all atoms within the radius
*
      If ((radiLD.gt.0.0d0).and.(nCtrLD.gt.0))  Then
#ifdef _MOLCAS_
         Call get_coord_all(Coord,nAtom)
#else
         call abend
#endif
         Do i=1,nCtrLD
           icen=iCtrLD(i)
           itmp1=3*icen-2
           Do j=1,nAtom
              If (j.ne.icen) Then
                 itmp2=3*j-2
                 distce=sqrt((Coord(itmp1)-Coord(itmp2))**2
     & +(Coord(itmp1+1)-Coord(itmp2+1))**2+(Coord(itmp1+2)-
     & Coord(itmp2+2))**2)
                 If (distce.le.radiLD) Then
                    indx2(j,3)=indx2(icen,3)
                 End If
              End If
           End Do
         End Do
*
**     Reorganize indx2
*
         Do i=1,nAtom-1
            Do j=i,nAtom
              If (indx2(i,3).gt.indx2(j,3)) Then
                Do k=1,3
                  itmp1=indx2(i,k)
                  indx2(i,k)=indx2(j,k)
                  indx2(j,k)=itmp1
                End Do
              End If
           End Do
         End Do
      End If
*
**    Compute the number of blocks and the maximum size
*
      Nsm=indx2(1,2)-indx2(1,1)+1
      nblock=1
      itmp1=Nsm
      indx2(1,4)=1
      Do j=2,nAtom
        If (indx2(j,3).ne.indx2(j-1,3)) Then
           nblock=nblock+1
           indx2(nblock,4)=j
           If (itmp1.gt.Nsm) Nsm=itmp1
           itmp1=0
        End If
        itmp1=itmp1+indx2(j,2)-indx2(j,1)+1
      End Do
      If (itmp1.gt.Nsm) Nsm=itmp1
      End
************************************************************************
************************************************************************
************************************************************************
      subroutine diag_ldkh(nAtom,nblock,N,isize,Nsm,indx2,ssm,hsm,s,h,
     *                     sinv,ewsm,eigsm,ew,eig,v,scr8,scrsq,scrt)
************************************************************************
*                                                                      *
*      Block-diagonalize the hamiltonian                               *
*                                                                      *
*                                                                      *
************************************************************************
      Implicit none
      Integer i,j,k,l,nblock,firstc,lastc,nAtom,indx2(nAtom,4),counter,
     *        jcen,icen,N,isize,Nsm,itrian
      Real*8 h(isize),s(isize),v(isize),ew(N),eig(N,N),
     *                 scr8(isize)
      Real*8 ssm(Nsm,Nsm),hsm(Nsm*(Nsm+1)/2),sinv(Nsm,Nsm),
     *                 scrt(Nsm*(Nsm+1)/2),ewsm(Nsm),eigsm(Nsm,Nsm),
     *                 scrsq(Nsm,Nsm)
*                                                                      *
************************************************************************
*                                                                      *
      Do i=1,isize
         scr8(i)=s(i)
      End Do
      Do i=1,N
        Do j=1,N
          eig(i,j)=0.0d0
        End Do
      End Do
*
      Do i=1,nblock
        firstc=indx2(i,4)
        if (i.ne.nblock) Then
           lastc=indx2(i+1,4)-1
        Else
           lastc=nAtom
        End If
        counter=0
        itrian=1
        Do jcen=firstc,lastc
          counter=counter+indx2(jcen,2)-indx2(jcen,1)+1
          Do k=indx2(jcen,1),indx2(jcen,2)
             Do icen=firstc,jcen-1
                Do j=indx2(icen,1),indx2(icen,2)
                  ssm(itrian,1)=scr8(k*(k-1)/2+j)
                  hsm(itrian)=h(k*(k-1)/2+j)
                  itrian=itrian+1
                End Do
             End Do
             Do j=indx2(jcen,1),k
               ssm(itrian,1)=scr8(k*(k-1)/2+j)
               hsm(itrian)=h(k*(k-1)/2+j)
               itrian=itrian+1
             End Do
          End Do
        End Do
*
**      Run the calculation with the small variables
*
        Call sog(counter,ssm,sinv,scrt,scrsq,ewsm)
        Call diag_dkh(hsm,counter,eigsm,ewsm,sinv,scrsq,1)
*
*      Reconstruct EIG and EW for the perturbation
*
        l=0
        itrian=1
        Do icen=firstc,lastc
           Do k=indx2(icen,1),indx2(icen,2)
              ew(k)=ewsm(l+1)
              Do j=indx2(icen,1),indx2(icen,2)
                eig(j,k)=eigsm(itrian,1)
                itrian=itrian+1
              End Do
              l=l+1
           End Do
        End Do
      End Do
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(v)
      End
************************************************************************
************************************************************************
************************************************************************
      subroutine ldkhpert(N,isize,h,s,sinv,scrsq,scrt,nblock,indx2,eig,
     *                    ew,nAtom,sprov)
************************************************************************
*                                                                      *
*     Compute the 2nd order perturbation for the Local DKH procedure   *
*     for eigenvectors and eigenvalues                                 *
*                                                                      *
*                                                                      *
************************************************************************
      use DKH_Info
      Implicit real*8(a-h,o-z)
      Integer N,isize,i,j,k,l,m,icen,nAtom,
     *        nblock,indx2(nAtom,4)
      Real*8 h(isize),s(isize),sinv(N,N),eig(N,N),ew(N)
      Real*8 scrsq(N,N,3),scrt(isize,2),norm
      Real*8 sprov(N,N)
*                                                                      *
************************************************************************
*                                                                      *
      Call mat_sq_from_t (scrsq(1,1,2),N,h)
      Call sog(N,s,sinv,scrt(1,1),scrt(1,2),scrsq(1,1,1))
      Call Mat_sq_from_t (sprov,N,s)
*
      Do i=1,N
         Do j=1,N
            scrsq(i,j,1)=0.0d0
            scrsq(i,j,3)=0.0d0
         End Do
         scrt(i,1)=ew(i)
      End Do
*
*          H*C
*
      Do i=1,nCtrLD
         icen=iCtrLD(i)
         Do j=1,nblock
         If (icen.ne.j) Then
            Do k=indx2(j,1),indx2(j,2)
               Do l=indx2(icen,1),indx2(icen,2)
                  scrsq(k,l,3)=0.0d0
                  scrsq(k,l,1)=0.0d0
                  Do m=indx2(icen,1),indx2(icen,2)
                     scrsq(k,l,3)=scrsq(k,l,3)+scrsq(k,m,2)*eig(m,l)
                     scrsq(k,l,1)=scrsq(k,l,1)+sprov(k,m)  *eig(m,l)
                  End Do
               End Do
            End Do
         End If
         End Do
      End Do
*
*          C*H
*
      Do i=1,nCtrLD
         icen=iCtrLD(i)
         Do j=1,nblock
         If (icen.ne.j) Then
            Do k=indx2(j,1),indx2(j,2)
               Do l=indx2(icen,1),indx2(icen,2)
                  scrsq(k,l,2)=0.0d0
                  sprov(k,l)  =0.0d0
                  Do m=indx2(j,1),indx2(j,2)
                     scrsq(k,l,2)=scrsq(k,l,2)+eig(m,k)*scrsq(m,l,3)
                     sprov(k,l)  =sprov(k,l)  +eig(m,k)*scrsq(m,l,1)
                  End Do
               End Do
            End Do
         End If
         End Do
      End Do
*
**    Perturbation
*
      Do i=1,nCtrLD
         icen=iCtrLD(i)
         Do j=1,nblock
         If (icen.ne.j) Then
            Do k=indx2(j,1),indx2(j,2)
               Do l=indx2(icen,1),indx2(icen,2)
                 denom=1d0/(ew(l)-ew(k))
                 norm=scrsq(k,l,2)-ew(l)*sprov(k,l)
                 scrt(l,1)=scrt(l,1)+(norm**2)*denom
                 eig(k,l)=norm*denom
               End Do
            End Do
         End If
         End Do
      End Do
*
**   Orthonormalization
*
      Do i=1,N
         ew(i)=scrt(i,1)
         Do j=i,N
            norm=0.0d0
            Do k=1,N
               norm=norm+eig(k,i)*eig(k,j)
            End Do
            if (i.eq.j) Then
              Do k=1,N
                 eig(k,i)=eig(k,i)/Sqrt(norm)
              End Do
            Else
              Do k=1,N
                 eig(k,j)=eig(k,j)-norm*eig(k,i)
              End Do
            End If
         End Do
      End Do
*
      End
