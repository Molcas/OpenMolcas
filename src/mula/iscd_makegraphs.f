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
      Subroutine ISCD_MakeGraphs(m_max,maxOrd,maxIncOrd,
     &                        Graph1,Graph2,nOsc)
C!
      Implicit Real*8 ( a-h,o-z )
      Integer Graph1 (m_max+1,nOsc+1)
      Integer Graph2 (m_max+1,m_max+1,nOsc)
      Integer nTabDim
#include "WrkSpc.fh"
C!
C!---- Initialize.
      If ( m_max.eq.0 ) Return
      Call TabDim_drv(m_max,nOsc,nTabDim)
      maxOrd = nTabDim-1
C!
C!---- Set up the vertex table
      do iv=1,m_max+1
        do jv=1,nOsc+1
          Graph1(iv,jv) = 0
        enddo
      enddo
      do iv=1,m_max+1
        Graph1(iv,2) = 1
      enddo
      do jv=1,nOsc+1
        Graph1(1,jv) = 1
      enddo
      If ( nOsc.gt.1 ) Then
      Do iOsc = 2,nOsc
        n = 0
        Do nQuanta = 0,m_max
          n = n+Graph1(nQuanta+1,iOsc)
           Graph1(nQuanta+1,iOsc+1) = n
        End Do
      End Do
      End If
C!
C!---- set up the arc table
      Call GetMem('Number','Allo','INTE',ipNumber,m_max+1)
      do iv=0,m_max
        iWork(ipNumber+iv) = 0
      enddo
      N = 0
      Do m = 1,m_max
        N = N+Graph1(m,nosc+1)
        iWork(ipNumber+m) = n
      End Do
      do iv=1,m_max+1
        do jv=1,m_max+1
          do kv=1,nOsc
            Graph2(iv,jv,kv) = 0
          enddo
        enddo
      enddo
      Do iOsc = 1,nosc
        Do iQ1 = 0,m_max         ! Where we are going
          Do iQ2 = 0,iQ1-1       ! Where we came from
            Do i = iQ2+1,iq1     ! Sum over preceding paths
              Graph2(iQ1+1,iQ2+1,iOsc) = Graph1(i+1,iOsc)+
     &             Graph2(iQ1+1,iQ2+1,iOsc)
            End Do
          End Do
        End Do
      End Do
C!
      Do iQ1 = 0,m_max            ! Where we are going
        Do iQ2 = 0,iq1            ! Where we came from
          Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+
     &       iWork(ipNumber+iQ1)
        End Do
      End Do
C!
      Call GetMem('Number','Free','INTE',ipNumber,m_max+1)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(maxIncOrd)
      End
