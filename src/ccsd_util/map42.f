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
       subroutine map42 (a,b,dimp,dimq,dimr,dims,dim1,dim2,dim3,dim4,p,
     &                   q,r,s,nfact)
c
       integer dimp,dimq,dimr,dims,dim1,dim2,dim3,dim4,p,q,r,s,nfact
       real*8 a(1:dimp,1:dimq,1:dimr,1:dims)
       real*8 b(1:dim1,1:dim2,1:dim3,1:dim4)
c     integer index(1:4)
c
       integer pp,qq,rr,ss
c
       if (nfact.eq.1) then
c
c     factor + 1
c
c     do 100 ss=1,dims
c     index(s)=ss
c     do 100 rr=1,dimr
c     index(r)=rr
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=a(pp,qq,rr,ss)
c100  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 111 ss=1,dims
       do 1110 rr=1,dimr
       do 1111 qq=1,dimq
       do 1112 pp=1,dimp
       b(pp,qq,rr,ss)=a(pp,qq,rr,ss)
 1112   continue
 1111   continue
 1110   continue
 111    continue
       else
c     124* (3)
       do 112 ss=1,dims
       do 1120 rr=1,dimr
       do 1121 qq=1,dimq
       do 1122 pp=1,dimp
       b(pp,qq,ss,rr)=a(pp,qq,rr,ss)
 1122   continue
 1121   continue
 1120   continue
 112    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 113 ss=1,dims
       do 1130 rr=1,dimr
       do 1131 qq=1,dimq
       do 1132 pp=1,dimp
       b(pp,rr,qq,ss)=a(pp,qq,rr,ss)
 1132   continue
 1131   continue
 1130   continue
 113    continue
       else
c     134* (2)
       do 114 ss=1,dims
       do 1140 rr=1,dimr
       do 1141 qq=1,dimq
       do 1142 pp=1,dimp
       b(pp,ss,qq,rr)=a(pp,qq,rr,ss)
 1142   continue
 1141   continue
 1140   continue
 114    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 115 ss=1,dims
       do 1150 rr=1,dimr
       do 1151 qq=1,dimq
       do 1152 pp=1,dimp
       b(pp,rr,ss,qq)=a(pp,qq,rr,ss)
 1152   continue
 1151   continue
 1150   continue
 115    continue
       else
c     143* (2)
       do 116 ss=1,dims
       do 1160 rr=1,dimr
       do 1161 qq=1,dimq
       do 1162 pp=1,dimp
       b(pp,ss,rr,qq)=a(pp,qq,rr,ss)
 1162   continue
 1161   continue
 1160   continue
 116    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 121 ss=1,dims
       do 1210 rr=1,dimr
       do 1211 qq=1,dimq
       do 1212 pp=1,dimp
       b(qq,pp,rr,ss)=a(pp,qq,rr,ss)
 1212   continue
 1211   continue
 1210   continue
 121    continue
       else
c     214* (3)
       do 122 ss=1,dims
       do 1220 rr=1,dimr
       do 1221 qq=1,dimq
       do 1222 pp=1,dimp
       b(qq,pp,ss,rr)=a(pp,qq,rr,ss)
 1222   continue
 1221   continue
 1220   continue
 122    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 123 ss=1,dims
       do 1230 rr=1,dimr
       do 1231 qq=1,dimq
       do 1232 pp=1,dimp
       b(rr,pp,qq,ss)=a(pp,qq,rr,ss)
 1232   continue
 1231   continue
 1230   continue
 123    continue
       else
c     234* (1)
       do 124 ss=1,dims
       do 1240 rr=1,dimr
       do 1241 qq=1,dimq
       do 1242 pp=1,dimp
       b(ss,pp,qq,rr)=a(pp,qq,rr,ss)
 1242   continue
 1241   continue
 1240   continue
 124    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 125 ss=1,dims
       do 1250 rr=1,dimr
       do 1251 qq=1,dimq
       do 1252 pp=1,dimp
       b(ss,pp,rr,qq)=a(pp,qq,rr,ss)
 1252   continue
 1251   continue
 1250   continue
 125    continue
       else
c     241* (3)
       do 126 ss=1,dims
       do 1260 rr=1,dimr
       do 1261 qq=1,dimq
       do 1262 pp=1,dimp
       b(rr,pp,ss,qq)=a(pp,qq,rr,ss)
 1262   continue
 1261   continue
 1260   continue
 126    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 131 ss=1,dims
       do 1310 rr=1,dimr
       do 1311 qq=1,dimq
       do 1312 pp=1,dimp
       b(qq,rr,pp,ss)=a(pp,qq,rr,ss)
 1312   continue
 1311   continue
 1310   continue
 131    continue
       else
c     314* (2)
       do 132 ss=1,dims
       do 1320 rr=1,dimr
       do 1321 qq=1,dimq
       do 1322 pp=1,dimp
       b(qq,ss,pp,rr)=a(pp,qq,rr,ss)
 1322   continue
 1321   continue
 1320   continue
 132    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 133 ss=1,dims
       do 1330 rr=1,dimr
       do 1331 qq=1,dimq
       do 1332 pp=1,dimp
       b(rr,qq,pp,ss)=a(pp,qq,rr,ss)
 1332   continue
 1331   continue
 1330   continue
 133    continue
       else
c     324* (1)
       do 134 ss=1,dims
       do 1340 rr=1,dimr
       do 1341 qq=1,dimq
       do 1342 pp=1,dimp
       b(ss,qq,pp,rr)=a(pp,qq,rr,ss)
 1342   continue
 1341   continue
 1340   continue
 134    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 135 ss=1,dims
       do 1350 rr=1,dimr
       do 1351 qq=1,dimq
       do 1352 pp=1,dimp
       b(rr,ss,pp,qq)=a(pp,qq,rr,ss)
 1352   continue
 1351   continue
 1350   continue
 135    continue
       else
c     342* (1)
       do 136 ss=1,dims
       do 1360 rr=1,dimr
       do 1361 qq=1,dimq
       do 1362 pp=1,dimp
       b(ss,rr,pp,qq)=a(pp,qq,rr,ss)
 1362   continue
 1361   continue
 1360   continue
 136    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 141 ss=1,dims
       do 1410 rr=1,dimr
       do 1411 qq=1,dimq
       do 1412 pp=1,dimp
       b(qq,ss,rr,pp)=a(pp,qq,rr,ss)
 1412   continue
 1411   continue
 1410   continue
 141    continue
       else
c     412* (3)
       do 142 ss=1,dims
       do 1420 rr=1,dimr
       do 1421 qq=1,dimq
       do 1422 pp=1,dimp
       b(qq,rr,ss,pp)=a(pp,qq,rr,ss)
 1422   continue
 1421   continue
 1420   continue
 142    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 143 ss=1,dims
       do 1430 rr=1,dimr
       do 1431 qq=1,dimq
       do 1432 pp=1,dimp
       b(rr,qq,ss,pp)=a(pp,qq,rr,ss)
 1432   continue
 1431   continue
 1430   continue
 143    continue
       else
c     423* (1)
       do 144 ss=1,dims
       do 1440 rr=1,dimr
       do 1441 qq=1,dimq
       do 1442 pp=1,dimp
       b(ss,qq,rr,pp)=a(pp,qq,rr,ss)
 1442   continue
 1441   continue
 1440   continue
 144    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 145 ss=1,dims
       do 1450 rr=1,dimr
       do 1451 qq=1,dimq
       do 1452 pp=1,dimp
       b(rr,ss,qq,pp)=a(pp,qq,rr,ss)
 1452   continue
 1451   continue
 1450   continue
 145    continue
       else
c     432* (1)
       do 146 ss=1,dims
       do 1460 rr=1,dimr
       do 1461 qq=1,dimq
       do 1462 pp=1,dimp
       b(ss,rr,qq,pp)=a(pp,qq,rr,ss)
 1462   continue
 1461   continue
 1460   continue
 146    continue
       end if
       end if
c
       end if
c
       else
c
c     factor = -1
c
c     do 200 ss=1,dims
c     index(s)=ss
c     do 200 rr=1,dimr
c     index(r)=rr
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=-a(pp,qq,rr,ss)
c200  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 211 ss=1,dims
       do 2110 rr=1,dimr
       do 2111 qq=1,dimq
       do 2112 pp=1,dimp
       b(pp,qq,rr,ss)=-a(pp,qq,rr,ss)
 2112   continue
 2111   continue
 2110   continue
 211    continue
       else
c     124* (3)
       do 212 ss=1,dims
       do 2120 rr=1,dimr
       do 2121 qq=1,dimq
       do 2122 pp=1,dimp
       b(pp,qq,ss,rr)=-a(pp,qq,rr,ss)
 2122   continue
 2121   continue
 2120   continue
 212    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 213 ss=1,dims
       do 2130 rr=1,dimr
       do 2131 qq=1,dimq
       do 2132 pp=1,dimp
       b(pp,rr,qq,ss)=-a(pp,qq,rr,ss)
 2132   continue
 2131   continue
 2130   continue
 213    continue
       else
c     134* (2)
       do 214 ss=1,dims
       do 2140 rr=1,dimr
       do 2141 qq=1,dimq
       do 2142 pp=1,dimp
       b(pp,ss,qq,rr)=-a(pp,qq,rr,ss)
 2142   continue
 2141   continue
 2140   continue
 214    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 215 ss=1,dims
       do 2150 rr=1,dimr
       do 2151 qq=1,dimq
       do 2152 pp=1,dimp
       b(pp,rr,ss,qq)=-a(pp,qq,rr,ss)
 2152   continue
 2151   continue
 2150   continue
 215    continue
       else
c     143* (2)
       do 216 ss=1,dims
       do 2160 rr=1,dimr
       do 2161 qq=1,dimq
       do 2162 pp=1,dimp
       b(pp,ss,rr,qq)=-a(pp,qq,rr,ss)
 2162   continue
 2161   continue
 2160   continue
 216    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 221 ss=1,dims
       do 2210 rr=1,dimr
       do 2211 qq=1,dimq
       do 2212 pp=1,dimp
       b(qq,pp,rr,ss)=-a(pp,qq,rr,ss)
 2212   continue
 2211   continue
 2210   continue
 221    continue
       else
c     214* (3)
       do 222 ss=1,dims
       do 2220 rr=1,dimr
       do 2221 qq=1,dimq
       do 2222 pp=1,dimp
       b(qq,pp,ss,rr)=-a(pp,qq,rr,ss)
 2222   continue
 2221   continue
 2220   continue
 222    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 223 ss=1,dims
       do 2230 rr=1,dimr
       do 2231 qq=1,dimq
       do 2232 pp=1,dimp
       b(rr,pp,qq,ss)=-a(pp,qq,rr,ss)
 2232   continue
 2231   continue
 2230   continue
 223    continue
       else
c     234* (1)
       do 224 ss=1,dims
       do 2240 rr=1,dimr
       do 2241 qq=1,dimq
       do 2242 pp=1,dimp
       b(ss,pp,qq,rr)=-a(pp,qq,rr,ss)
 2242   continue
 2241   continue
 2240   continue
 224    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 225 ss=1,dims
       do 2250 rr=1,dimr
       do 2251 qq=1,dimq
       do 2252 pp=1,dimp
       b(ss,pp,rr,qq)=-a(pp,qq,rr,ss)
 2252   continue
 2251   continue
 2250   continue
 225    continue
       else
c     241* (3)
       do 226 ss=1,dims
       do 2260 rr=1,dimr
       do 2261 qq=1,dimq
       do 2262 pp=1,dimp
       b(rr,pp,ss,qq)=-a(pp,qq,rr,ss)
 2262   continue
 2261   continue
 2260   continue
 226    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 231 ss=1,dims
       do 2310 rr=1,dimr
       do 2311 qq=1,dimq
       do 2312 pp=1,dimp
       b(qq,rr,pp,ss)=-a(pp,qq,rr,ss)
 2312   continue
 2311   continue
 2310   continue
 231    continue
       else
c     314* (2)
       do 232 ss=1,dims
       do 2320 rr=1,dimr
       do 2321 qq=1,dimq
       do 2322 pp=1,dimp
       b(qq,ss,pp,rr)=-a(pp,qq,rr,ss)
 2322   continue
 2321   continue
 2320   continue
 232    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 233 ss=1,dims
       do 2330 rr=1,dimr
       do 2331 qq=1,dimq
       do 2332 pp=1,dimp
       b(rr,qq,pp,ss)=-a(pp,qq,rr,ss)
 2332   continue
 2331   continue
 2330   continue
 233    continue
       else
c     324* (1)
       do 234 ss=1,dims
       do 2340 rr=1,dimr
       do 2341 qq=1,dimq
       do 2342 pp=1,dimp
       b(ss,qq,pp,rr)=-a(pp,qq,rr,ss)
 2342   continue
 2341   continue
 2340   continue
 234    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 235 ss=1,dims
       do 2350 rr=1,dimr
       do 2351 qq=1,dimq
       do 2352 pp=1,dimp
       b(rr,ss,pp,qq)=-a(pp,qq,rr,ss)
 2352   continue
 2351   continue
 2350   continue
 235    continue
       else
c     342* (1)
       do 236 ss=1,dims
       do 2360 rr=1,dimr
       do 2361 qq=1,dimq
       do 2362 pp=1,dimp
       b(ss,rr,pp,qq)=-a(pp,qq,rr,ss)
 2362   continue
 2361   continue
 2360   continue
 236    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 241 ss=1,dims
       do 2410 rr=1,dimr
       do 2411 qq=1,dimq
       do 2412 pp=1,dimp
       b(qq,ss,rr,pp)=-a(pp,qq,rr,ss)
 2412   continue
 2411   continue
 2410   continue
 241    continue
       else
c     412* (3)
       do 242 ss=1,dims
       do 2420 rr=1,dimr
       do 2421 qq=1,dimq
       do 2422 pp=1,dimp
       b(qq,rr,ss,pp)=-a(pp,qq,rr,ss)
 2422   continue
 2421   continue
 2420   continue
 242    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 243 ss=1,dims
       do 2430 rr=1,dimr
       do 2431 qq=1,dimq
       do 2432 pp=1,dimp
       b(rr,qq,ss,pp)=-a(pp,qq,rr,ss)
 2432   continue
 2431   continue
 2430   continue
 243    continue
       else
c     423* (1)
       do 244 ss=1,dims
       do 2440 rr=1,dimr
       do 2441 qq=1,dimq
       do 2442 pp=1,dimp
       b(ss,qq,rr,pp)=-a(pp,qq,rr,ss)
 2442   continue
 2441   continue
 2440   continue
 244    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 245 ss=1,dims
       do 2450 rr=1,dimr
       do 2451 qq=1,dimq
       do 2452 pp=1,dimp
       b(rr,ss,qq,pp)=-a(pp,qq,rr,ss)
 2452   continue
 2451   continue
 2450   continue
 245    continue
       else
c     432* (1)
       do 246 ss=1,dims
       do 2460 rr=1,dimr
       do 2461 qq=1,dimq
       do 2462 pp=1,dimp
       b(ss,rr,qq,pp)=-a(pp,qq,rr,ss)
 2462   continue
 2461   continue
 2460   continue
 246    continue
       end if
       end if
c
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(s)
       end
