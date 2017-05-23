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
c*****************************************************************
      subroutine complete_ext_loop()
c*****************************************************************
c     26 feb 2007 - revised
c
#include "drt_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei

c     write(6,*) ' ext_test '
c      return
      do irot=1,mcroot
        irtidx=indx(irot)
        lwei=isegsta !+irtidx

        do iupwei=1,isegupwei
          ilpvalue=0
          mm0=lwei
          nn=lwei+icano_nnsta-1
          do nntmp=icano_nnsta,icano_nnend
            nn=nn+1
            mm=mm0
            vlptmp1=vector1(nn+irtidx)
            vlptmp=vector2(nn+irtidx)
            do mmtmp=1,nntmp-1
              ilpvalue=ilpvalue+1
              vetmp=value_lpext(ilpvalue)
              mm=mm+1
              vector2(mm+irtidx)=vector2(mm+irtidx)+vlptmp1*vetmp
              vlptmp=vlptmp+vector1(mm+irtidx)*vetmp
            enddo
            vector2(nn+irtidx)=vlptmp
          enddo
          lwei=lwei+isegdownwei
        enddo

      enddo
c...end of complete_ext_loop
      end

      subroutine complete_ext_loop_g()
#include "drt_h.fh"
#include "grad_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)

      lwei=isegsta
      do iupwei=1,isegupwei
        ilpvalue=0
        mm0=lwei
        nn=lwei+icano_nnsta-1
        do nntmp=icano_nnsta,icano_nnend
         nn=nn+1
         mm=mm0
         do mmtmp=1,nntmp-1
            ilpvalue=ilpvalue+1
            mm=mm+1
            indexlp=index_lpext(ilpvalue)
          if(indexlp.ne.0) then
            valuelp=value_lpext(ilpvalue)
            vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
          end if
            indexlp1=index_lpext1(ilpvalue)
          if(indexlp1.ne.0) then
            valuelp1=value_lpext1(ilpvalue)
            vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
          end if
            indexlp2=index_lpext2(ilpvalue)
          if(indexlp2.ne.0) then
            valuelp2=value_lpext2(ilpvalue)
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
          endif
         enddo
        enddo
        lwei=lwei+isegdownwei
      enddo
      end
