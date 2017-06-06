PROGRAM cylindrical
implicit none
integer::j,k,l,m,n
integer,parameter::kmax=1000,lmax=1000
complex(kind(0d0))::T(1:3,1:3),det,detin,detinold
complex(kind(0d0)),parameter::i=dcmplx(0d0,1d0)
real(8)::C11,C12,C44,k0,k1,k2,q,qmax,Vt,Vl,omega,omegamax,rho,R,diff,bessjn,bessin,bessj0,bessj1
external bessjn,bessin,diff,bessj0,bessj1
OPEN(10,file='cyl2.dat')
!===================Cij(弾性定数),rho(密度),R(半径),Vl,Vt(速度)================
n=2
C11=3d11!11.88d11
C44=1d11!5.94d11
C12=C11-2d0*C44
rho=1d0!5.36d0
R=1d-6
Vl=DSQRT(C11/rho)
Vt=DSQRT(C44/rho)
omegamax=1d13
!===============無次元化====================
C11=C11/C44
C12=C12/C44
C44=C44/C44
rho=rho/rho
omegamax=omegamax*R/Vl
Vt=Vt/Vl
Vl=Vl/Vl
qmax=omegamax/Vl
R=R/R
!==================================================
do l=1,lmax
q=dble(l)/dble(lmax)*qmax
detin=0d0

 do k=1,kmax
   omega=dble(k)/dble(kmax+1)*(omegamax-q*Vl)+q*Vl
   k0=DSQRT((omega/Vl)**2-q**2)
   k1=DSQRT((omega/Vt)**2-q**2)
   k2=k1
!==============================================================================
T(1,1)=(-1d0/2d0*omega**2-k0**2+dble(n**2-n))*bessjn(n,k0)+k0*bessjn(n+1,k0)
T(1,2)=i*dble(n*(n-1))*bessjn(n,k1)-i*dble(n)*k1*bessjn(n+1,k1)
T(1,3)=-i*q*(k2**2+dble(n-n**2))/k2*bessjn(n,k2)+i*q*bessjn(n+1,k2)
T(2,1)=2d0*i*dble(n*(n-1))*bessjn(n,k0)-2d0*i*dble(n)*k0*bessjn(n+1,k0)
T(2,2)=(k1**2+2d0*dble(-n**2+n))*bessjn(n,k1)-2d0*k1*bessjn(n+1,k1)
T(2,3)=-2d0*q/k2*dble(n*(n-1))*bessjn(n,k2)+2d0*dble(n)*q*bessjn(n+1,k2)
T(3,1)=-2d0*i*q*k0*bessjn(n+1,k0)+2d0*i*q*dble(n)*bessjn(n,k0)
T(3,2)=-q*dble(n)*bessjn(n,k1)
T(3,3)=(2d0*q**2-omega**2*(VL/Vt)**2)*(bessjn(n+1,k2)-dble(n)/k2*bessjn(n,k2))
!==============================================================================
 call detN(3,T,det) 
 if(dble(det)*dble(detin)<0d0) then
 write(10,*) q,omega
end if
detin=1d0/det
end do
detin=0d0

 do k=1,kmax
   omega=dble(k)/dble(kmax+1)*q*(Vl-Vt)+q*Vt
   k0=DSQRT(q**2-(omega/Vl)**2)
   k1=DSQRT((omega/Vt)**2-q**2)
   k2=k1
!==============================================================================
T(1,1)=(-1d0/2d0*omega**2+k0**2+dble(n**2-n))*bessin(n,k0)-k0*bessin(n+1,k0)
T(1,2)=i*dble(n*(n-1))*bessjn(n,k1)-i*dble(n)*k1*bessjn(n+1,k1)
T(1,3)=-i*q*(k2**2+dble(n-n**2))/k2*bessjn(n,k2)+i*q*bessjn(n+1,k2)
T(2,1)=2d0*i*dble(n*(n-1))*bessin(n,k0)+2d0*i*dble(n)*k0*bessin(n+1,k0)
T(2,2)=(k1**2+2d0*dble(-n**2+n))*bessjn(n,k1)-2d0*k1*bessjn(n+1,k1)
T(2,3)=-2d0*q/k2*dble(n*(n-1))*bessjn(n,k2)+2d0*dble(n)*q*bessjn(n+1,k2)
T(3,1)=2d0*i*q*k0*bessin(n+1,k0)+2d0*i*dble(n)*q*bessin(n,k0)
T(3,2)=-q*dble(n)*bessjn(n,k1)
T(3,3)=(2d0*q**2-omega**2*(VL/Vt)**2)*(bessjn(n+1,k2)-dble(n)/k2*bessjn(n,k2))
!==============================================================================
  call detN(3,T,det)
 if(dble(det)*dble(detin)<0d0) then

 write(10,*) q,omega
end if
detin=1d0/det
end do
detin=0d0

do k=1,kmax
omega=dble(k)/dble(kmax+1)*q*Vt
   k0=DSQRT(q**2-(omega/Vl)**2)
   k1=DSQRT(q**2-(omega/Vt)**2)
   k2=k1
!==============================================================================
T(1,1)=(-1d0/2d0*omega**2+k0**2+dble(n**2-n))*bessin(n,k0)-k0*bessin(n+1,k0)
T(1,2)=i*dble(n*(n-1))*bessin(n,k1)+i*dble(n)*k1*bessin(n+1,k1)
T(1,3)=i*q*(k2**2+dble(n**2-n))/k2*bessin(n,k2)-i*q*bessin(n+1,k2)
T(2,1)=2d0*i*dble(n*(n-1))*bessin(n,k0)+2d0*i*dble(n)*k0*bessin(n+1,k0)
T(2,2)=(-k1**2+2d0*dble(-n**2+n))*bessin(n,k1)+2d0*k1*bessin(n+1,k1)
T(2,3)=-2d0*q/k2*dble(n*(n-1))*bessin(n,k2)-2d0*q*dble(n)*bessin(n+1,k2)
T(3,1)=2d0*i*q*k0*bessin(n+1,k0)+2d0*i*dble(n)*q*bessin(n,k0)
T(3,2)=-q*dble(n)*bessin(n,k1)
T(3,3)=-(2d0*q**2-omega**2*(VL/Vt)**2)*(dble(n)/k2*bessin(n,k2)+bessin(n+1,k2))
!==============================================================================
call detN(3,T,det)
if(dble(det)*dble(detin)<0d0) then
 write(10,*) q,omega
end if
detin=1d0/det
end do
detin=0d0
end do
END

!==========行列の固有値を計算する==============================================

subroutine detN(N,PO,w)
IMPLICIT NONE
integer::i,j,k,N
complex*16 Q(N,N),PO(N,N),P(N,N),w,AR

AR=(1d0,0d0)

do i=1,n
do j=1,n
P(i,J)=PO(i,j)
end do
end do

do k=1,N-1
do i=1,N
do j=1,N
Q(i,j) = P(i,j)
end do
end do
do i= k+1,N
do j=k,N
P(i,j)=P(i,j)-P(k,j)*Q(i,k)/Q(k,k)
end do
end do
end do
w=AR
do i=1,N
w=w*P(i,i)
end do
return
end subroutine

!==============整数n次のベッセル関数、変形ベッセル関数＝＝＝＝＝＝＝＝

double precision FUNCTION bessjn(n,x)
integer::n
real(8)::x,bessj0,bessj1,bessj
external bessj0,bessj1,bessj
 if (n.ge.2) then
   bessjn=bessj(n,x)
  else if (n==1) then
   bessjn=bessj1(x)
  else if (n==0) then
   bessjn=bessj0(x)
  else if (n==-1) then
   bessjn=-bessj1(x)
  else if (n.le.-2) then
   bessjn=(-1d0)**abs(n)*bessj(abs(n),x)
 end if
end function


double precision FUNCTION bessin(n,x)
integer::n
real(8)::x,bessi0,bessi1,bessi
external bessi0,bessi1,bessi
 if(n==0) then
   bessin=bessi0(x)
  else if (n==1.or.n==-1) then
   bessin=bessi1(x)
  else if (n.gt.1.or.n.lt.-1) then
   bessin=bessi(abs(n),x)
 end if
end function

!==============(n>1)のベッセル関数、変形ベッセル関数==================

      FUNCTION bessj(n,x)
      INTEGER::n,IACC
      REAL(8)::bessj,x,BIGNO,BIGNI
      PARAMETER(IACC=40,BIGNO=1.d10,BIGNI=1.d-10)
      external bessj0,bessj1
      INTEGER:: j,jsum,m
      REAL(8):: ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
      if(n.lt.2)pause 'bad argument n in bessj'
      ax=abs(x)
      if(ax.eq.0.d0)then
        bessj=0.
      else if(ax.gt.dble(n))then
        tox=2.d0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2.0d0/ax
        m=2*((n+int(dsqrt(dble(IACC*n))))/2)
        bessj=0.0d0
        jsum=0
        sum=0.0d0
        bjp=0.0d0
        bj=1.0d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
      return
      END


      FUNCTION bessi(n,x)
      INTEGER::n,IACC
      REAL(8)::bessi,x,BIGNO,BIGNI
      PARAMETER(IACC=40,BIGNO=1.0d10,BIGNI=1.0d-10)
      external bessi0
      INTEGER:: j,m
      REAL(8):: bi,bim,bip,tox,bessi0
      if (n.lt.2) pause 'bad argument n in bessi'
      if (x.eq.0.) then
        bessi=0.
      else
        tox=2.0/abs(x)
        bip=0.0
        bi=1.0
        bessi=0.
        m=2*((n+int(dsqrt(dble(IACC*n)))))
        do 11 j=m,1,-1
          bim=bip+dble(j)*tox*bi
          bip=bi
          bi=bim
          if (abs(bi).gt.BIGNO) then
            bessi=bessi*BIGNI
            bi=bi*BIGNI
            bip=bip*BIGNI
          endif
          if (j.eq.n) bessi=bip
11      continue
        bessi=bessi*bessi0(x)/bi
        if (x.lt.0..and.mod(n,2).eq.1) bessi=-bessi
      endif
      return
      END

!===================0次,1次のベッセル、変形ベッセル============================

      FUNCTION bessj0(x)
      REAL(8):: bessj0,x
      REAL(8):: ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.0d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=dsqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END


      FUNCTION bessj1(x)
      REAL(8):: bessj1,x
      REAL(8):: ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.0d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=dsqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      END


      FUNCTION bessi0(x)
      REAL(8):: bessi0,x
      REAL(8):: ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,-0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75d0)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi0=(exp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      END


      FUNCTION bessi1(x)
      REAL(8):: bessi1,x
      REAL(8):: ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75d0)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi1=(exp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        if(x.lt.0.)bessi1=-bessi1
      endif
      return
      END
