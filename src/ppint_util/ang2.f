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
      subroutine ang2_molcas(ang,binom,crda,dfac,it,l,lit,
     &  lmlo,lmhi,lmf,lml,
     &  lmx,lmy,lmz,lmnv,mproju,xk,yk,zk,zlm)
c
c  compute type 2 angular integrals
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0)
      dimension ang(lit,mproju,*), binom(*), crda(lit,3), dfac(*),
     &  lmf(*), lml(*), lmx(*), lmy(*), lmz(*), lmnv(3,*), zlm(*)
c
      call wzero(lit*mproju*lmhi,ang,1)
      na1=lmnv(1,it)+1
      la1=lmnv(2,it)+1
      ma1=lmnv(3,it)+1
      naind=(na1*(na1-1))/2
      laind=(la1*(la1-1))/2
      maind=(ma1*(ma1-1))/2
      loc1=(l-1)**2
      mhi=l+l-1
      do 80 ia=1,na1
      pab1=binom(naind+ia)*crda((na1+1)-ia,1)
      if(pab1.eq.a0) go to 80
      do 70 ib=1,la1
      pab2=pab1*binom(laind+ib)*crda((la1+1)-ib,2)
      if(pab2.eq.a0) go to 70
      do 60 ic=1,ma1
      pab3=pab2*binom(maind+ic)*crda((ma1+1)-ic,3)
      if(pab3.eq.a0) go to 60
      n=((ia-3)+ib)+ic
      lamlo=max(l-n,lmlo+mod(l+n+lmlo,2))
      lamhi=min(l+n,lmhi-mod(l+n+lmhi,2))
      if(lamlo.gt.lamhi) go to 60
      do 55 m=1,mhi
      mstart=lmf(loc1+m)
      mend=lml(loc1+m)
      do 50 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc2=(lam-1)**2
      do 40 mu=1,l2
      istart=lmf(loc2+mu)
      if(mod(ia+lmx(mstart)+lmx(istart),2).ne.1.or.
     1   mod(ib+lmy(mstart)+lmy(istart),2).ne.1.or.
     2   mod(ic+lmz(mstart)+lmz(istart),2).ne.1) go to 40
      pre=a0
      iend=lml(loc2+mu)
      aint=a0
      do 30 i=istart,iend
        indx=lmx(i)
        indy=lmy(i)
        indz=lmz(i)
        if(indx.eq.0) then
          xkp=a1
        else
          xkp=xk**indx
        endif
        if(indy.eq.0) then
          ykp=a1
        else
          ykp=yk**indy
        endif
        if(indz.eq.0) then
          zkp=a1
        else
          zkp=zk**indz
        endif
        pre=pre+zlm(i)*xkp*ykp*zkp
        do 20 j=mstart,mend
          mndx=lmx(j)
          mndy=lmy(j)
          mndz=lmz(j)
          aint=aint+zlm(i)*zlm(j)*
     &              dfac(ia+indx+mndx)*
     &              dfac(ib+indy+mndy)*
     &              dfac(ic+indz+mndz)/
     &              dfac(ia+indx+mndx+ib+indy+mndy+ic+indz+mndz)
   20   continue
   30 continue
      angt=angt+pre*aint
   40 continue
      ang(n+1,m,lam)=ang(n+1,m,lam)+angt*pab3
   50 continue
   55 continue
   60 continue
   70 continue
   80 continue
      return
      end
