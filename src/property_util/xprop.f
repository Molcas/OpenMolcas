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
      Subroutine Xprop(short,
     &                 nirrep,nbas,ntotv,vec,ntoto,occ,thrs,
     &                 ntotd,opel,
     &                 out)
************************************************************************
*                                                                      *
*     Purpose: the calculation of the average value of an operator,    *
*              whose matrix elements are supplied in the array opel,   *
*              for a set of eigenvectors supplied in vec with occupa-  *
*              tion numbers supplied in occ                            *
*                                                                      *
*     Parameters                                                       *
*                                                                      *
*       short          logical, if (short) than only the total average *
*                      value is calculated and transferred in out(1)   *
*                                                                      *
*       nirrep         number of irreps                                *
*       nbas(0:nirrep) dimension for each irrep                        *
*                                                                      *
*       ntotv          the total number of elemnts for all eigenvectors*
*       vec            if (short) then vec stores all lower triangles  *
*       (1:ntotv)      for all diagonal blocks of the density matrix   *
*                      size: sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2)*
*                      else                                            *
*                      vec stores the eigenvectors                     *
*                      size: sum(i,i=0,nirrep-1)(nbas(i)*nbas(i))      *
*       ntoto          the total number of vectors for all represen-   *
*                      tations=the number of basis functions           *
*       occ            if (short) the occ array is a dummy             *
*       (1:ntoto)      else                                            *
*                      occ stores the occupation numbers               *
*                      size: sum(i,i=0,nirrep-1)(nbas(i))              *
*       thrs           if (short) then this parameter is a dummy       *
*                      else                                            *
*                      if the orbital occupation number is .le.        *
*                      thrs the orbital contribution will not be       *
*                      printed out.                                    *
*                                                                      *
*       ntotd          the total number of elements in lower triangles *
*                      of all diagonal blocks                          *
*       opel           a storage area for transferring all lower       *
*       (1:ntotd)      triangles of diagonal blocks of the operator    *
*                      matrix                                          *
*                      size: sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2)*
*       out            on return if (short) out(1) contains the total  *
*                      average value                                   *
*                      else                                            *
*                      out(i), i=1,sum(k,k=0,nirrep-1)(nbas(i))        *
*                      contains the orbital contributions (multiplied  *
*                      by the corresponding occupation numbers)        *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      logical short
      dimension nbas(0:nirrep-1),vec(1:ntotv),
     &          occ(1:ntoto),opel(1:ntotd),out(1:ntoto)
*
      if (short) then
        icount=1
        sum=DDOT_(ntotd,vec,1,opel,1)
        Out(1) = Sum
      else
        ndim2=0
        do 1 i=0,nirrep-1
          ndim2=ndim2+nbas(i)**2
    1   continue
*
        iadv=0
        iado=0
        iadout=0
        jCount = 1
        do 199 i=0,nirrep-1
        do 200 iv=1,nbas(i)
          iado=iado+1
          iadout=iadout+1
          sum=0.0d+00
          icount=jCount
          do 201 iv1=1,nbas(i)
            do 202 iv2=1,iv1-1
              sum=sum+2.0d+00*vec(iadv+iv1)*vec(iadv+iv2)*opel(icount)
              icount=icount+1
  202       continue
            sum=sum+vec(iadv+iv1)*vec(iadv+iv1)*opel(icount)
            icount=icount+1
  201     continue
          out(iadout)=occ(iado)*sum
          iadv=iadv+nbas(i)
  200   continue
        jCount = jCount + nBas(i)*(nBas(i)+1)/2
  199   continue
      endif
*
*     if (short) then
*       write (*,'(3x,a,f18.10)') 'Total = ', out(1)
*     else
*       ii=0
*       do 2 i=0,nirrep-1
*         write (*,'(1x,a,i2)') 'Irrep No.',i
*         write (*,'(5(3x,i3,i3,f18.10))')
*    &    (j,out(ii+j),j=1,nbas(i))
*         ii=ii+nbas(i)
*   2   continue
*     endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(thrs)
      End
