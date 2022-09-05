!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Mk_List2(List2,nTheta_All,mData,nSO_Tot,iCnttp,nTest,  &
     &                    ijS_req)
      Use Basis_Info, only: dbsc, Shells
#include "stdalloc.fh"
      Integer List2(2*mData,nTheta_All)
      Logical Only_DB

      Integer, allocatable:: iList(:,:)
!
      Call mma_allocate(iList,mData,nSO_Tot,Label='iList')
!
      Only_DB=ijS_req.ne.0
!
!     Generate intermediate list
!     Generate list2, shell blocked!
!
      ijSO=0
      iiSO=0
      iSO_= 0
      Do iAng = 0, nTest
         iShll = dbsc(iCnttp)%iVal + iAng
         nCmp = (iAng+1)*(iAng+2)/2
         If (Shells(iShll)%Prjct) nCmp = 2*iAng+1
         nSO=nCmp*Shells(iShll)%nBasis
         Do iCmp = 1, nCmp
            nCont = Shells(iShll)%nBasis
            Do iCont = 1, nCont
                iSO_= iSO_+ 1
                iList(1,iSO_)=iAng
                iList(2,iSO_)=iCmp
                iList(3,iSO_)=iCont
                iList(4,iSO_)=iShll
            End Do
         End Do
!        Write (6,*) 'iSO_=',iSO_
!
         jjSO=0
         Do jAng = 0, iAng
!           Write (6,*) 'iAng,jAng=',iAng,jAng
            jShll = dbsc(iCnttp)%iVal + jAng
            mCmp = (jAng+1)*(jAng+2)/2
            If (Shells(jShll)%Prjct) mCmp = 2*jAng+1
!
            mSO=mCmp*Shells(jShll)%nBasis
!
            ijS=(iAng+1)*iAng/2+jAng+1
!
            If (.NOT.Only_DB .or. ijS.eq.ijS_req) Then
               Do iSO = iiSO+1, iiSO+nSO
                  iAng_ =iList(1,iSO)
                  iCmp_ =iList(2,iSO)
                  iCont_=iList(3,iSO)
                  iShll_=iList(4,iSO)
!
                  jSO_Max=jjSO+mSO
                  If (jAng.eq.iAng) jSO_Max=iSO
                  Do jSO = jjSO+1, jSO_Max
                     ijSO=ijSO+1
                     jAng_ =iList(1,jSO)
                     jCmp_ =iList(2,jSO)
                     jCont_=iList(3,jSO)
                     jShll_=iList(4,jSO)
!
!                    Write (*,*) 'iSO,jSO,ijSO=',iSO,jSO,ijSO
                     List2(1,ijSO)=iAng_
                     List2(2,ijSO)=jAng_
                     List2(3,ijSO)=iCmp_
                     List2(4,ijSO)=jCmp_
                     List2(5,ijSO)=iCont_
                     List2(6,ijSO)=jCont_
                     List2(7,ijSO)=iShll_
                     List2(8,ijSO)=jShll_
                  End Do            ! jSO
               End Do               ! iSO
!
            End If
!
            jjSO=jjSO+mSO
         End Do                  ! jAng
!
         iiSO=iiSO+nSO
      End Do                     ! iAng
!
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
         Write (6,*) 'List2'
         Write (6,*) '  iAng,  jAng,  iCmp,  jCmp, iCont, '             &
     &             //'jCont, iShll, jShll'
         Do ijSO = 1, nTheta_All
            Write (6,'(8I7)') (List2(i,ijSO),i=1,8)
         End Do
#endif
!
      Call mma_deallocate(iList)
!
      Return
      End
