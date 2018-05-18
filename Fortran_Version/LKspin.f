!     F. Antonini and C. Rodriguez; May, 17, 2018 
!     Variable step-size Runge-Kutta-Fehlberg method 7(8) is used
      
      module my_subs
 
      implicit none

      contains
            
      FUNCTION cross(a, b)
      double precision, DIMENSION(3) :: cross
      double precision, DIMENSION(3), INTENT(IN) :: a, b
      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
      END FUNCTION cross

      FUNCTION dotp(a, b)
      double precision dotp
      double precision, DIMENSION(3), INTENT(IN) :: a, b
      dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      END FUNCTION dotp


      SUBROUTINE random(rnd)
!     genera un numero casuale tra 0 ed 1
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: rnd
      INTEGER :: num = 10
      INTEGER :: isize,idate(8)
      INTEGER,ALLOCATABLE :: iseed(:)
      INTEGER :: i
      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize) !imposta il numero di interi per contenere il seme      
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed) !acquisisce il valore corrente del seme
      iseed = iseed * (idate(8)-1000) ! idate(8) contains milisecond
      CALL RANDOM_SEED(PUT=iseed) !imposta il nuovo seme            
      CALL RANDOM_NUMBER(rnd)
      DEALLOCATE(iseed)           
      END SUBROUTINE random

      
********************************************************************
C     CALCULATE TIME DERIVATIVES                                    C
********************************************************************
      subroutine force(t,y,yp)
      implicit double precision (a-z)
      integer n,k,h
      parameter  (n=19)
      double precision y(n),yp(n),S0(3),Seff(3)
      double precision j(3),edot(3),jdot(3),e2dot(3),j2dot(3)
      double precision S1dot(3),S2dot(3),S1(3),S2(3)
      double precision ev(3),jv(3),ev2(3),jv2(3),lai(3)
      common/cnst/c,pi
      common/orbit/m1,m2,m3,a2,a1
      common/par/GR,SP,SS,OC,OR
      yp=0.
      
      do h=1,3
         ev(h)=y(h)
      end do
      do h=1,3
         jv(h)=y(h+3)
      end do
      do h=1,3
         ev2(h)=y(h+6)
      end do
      do h=1,3
         jv2(h)=y(h+9)
      end do
      
!     define various quantities      
      a1=y(19)
      e1=sqrt(y(1)**2+y(2)**2+y(3)**2)
      e12=(y(1)**2+y(2)**2+y(3)**2)    
      j1=sqrt(y(4)**2+y(5)**2+y(6)**2)
      j2t=sqrt(y(10)**2+y(11)**2+y(12)**2)
      e2=sqrt(y(7)**2+y(8)**2+y(9)**2)
      
      P1=2.*pi*sqrt(a1**3/(m1+m2))
      mm=sqrt((m1+m2)/a1**3)
      
      mu1=m1*m2/(m1+m2)
      mu2=(m1+m2)*m3/(m1+m2+m3)
      L1=mu1*sqrt((m1+m2)*a1)
      L2=mu2*sqrt((m1+m2+m3)*a2)
         
!     KOZAI ts
      TLK=(m1+m2)*(a2*j2t/a1)**3/mm/m3*4./3.

!     GR terms      
      om_gr=3.*mm*(m1+m2)/c**2/(a1*j1**2)*GR
      e_GW=-304.*m1*m2*(m1+m2)*e1*(1.+e1**2*121./304.)
     &     /(15.*c**5*a1**4*((1.-e1**2)**(5./2.)))*GR
      l_GW=-e_GW*e1/sqrt(1.-e1**2)*GR

!     Some dot-prod
      e1j2=dotp(ev,jv2/j2t)
      j1j2=dotp(jv,jv2/j2t)
      e1e2=dotp(ev,ev2/e2)
      j1e2=dotp(jv,ev2/e2)
     
!     SMA evolution
      yp(19)=(-64.*m1*m2*(m1+m2)/(5.*(c**5)*(y(19)**3)*(1-e1**2)**3.5))
     &     *(1.+(e1**2)*73./24.+(e1**4)*37./96.)*GR
      
!     EV of INNER ORBIT
      oct=(m1-m2)/(m1+m2)*a1/a2*e2/j2t**2*75./64.*(4./3.)*OC
           
      edot=(2.d0*cross(jv,ev)
     &     -5.d0*e1j2*cross(jv,jv2/j2t)
     &     +j1j2*cross(ev,jv2/j2t))/TLK 
     &     +om_gr*cross(jv/j1,ev) !GR precession
     &     +e_GW*ev/e1          !GW
     &     -(2.*e1j2*j1j2*cross(ev,ev2/e2)
     &     +(8./5.*e1**2-1./5.-7.*e1j2**2+j1j2**2)*cross(jv,ev2/e2)
     &     +2.*(e1e2*j1j2+e1j2*j1e2)*cross(ev,jv2/j2t)
     &     +2.*(j1j2*j1e2-7.*e1j2*e1e2)*cross(jv,jv2/j2t)
     &     +16./5.*e1e2*cross(jv,ev))
     &     *oct/TLK
      
      jdot=-(2.*(e1e2*j1j2+e1j2*j1e2)*cross(jv,jv2/j2t)
     &     +2.*(j1e2*j1j2-7.*e1e2*e1j2)*cross(ev,jv2/j2t)
     &     +2.*(e1j2*j1j2)*cross(jv,ev2/e2)
     &     +(8./5.*e1**2-1./5.-7.*e1j2**2+j1j2**2)*cross(ev,ev2/e2))
     &     *oct/TLK
     &     +(j1j2*cross(jv,jv2/j2t)
     &     -5.d0*e1j2*cross(ev,jv2/j2t))/TLK !djx/dt
     &     +l_GW*jv/j1          !GW

      yp(1)=edot(1)
      yp(2)=edot(2)
      yp(3)=edot(3)
       
      yp(4)=jdot(1)      
      yp(5)=jdot(2)
      yp(6)=jdot(3)


!     EVOLUTION OF OUTER ORBIT  / OCTUPOLE
      e2dot=(j1j2*cross(ev2,jv)
     &     -5.*e1j2*cross(ev2,ev)
     &     -(1./2.-3.*e1**2+25./2.*e1j2**2
     &     -5./2.*j1j2**2)*cross(jv2/j2t,ev2))
     &     /TLK*L1/L2/j2t
     &     -(2.*(e1j2*j1e2*cross(ev2/e2,jv)
     &     +j1j2*e1e2*cross(ev2/e2,jv)
     &     +j2t**2/e2*e1j2*j1j2*cross(jv2/j2t,jv))
     &     +2.*j1e2*j1j2*cross(ev2/e2,ev)
     &     -14.*e1e2*e1j2*cross(ev2/e2,ev)
     &     +j2t**2/e2*(8./5.*e1**2-1./5.
     &     -7.*e1j2**2+j1j2**2)*cross(jv2/j2t,ev)
     &     -2.*(1./5.-8./5.*e1**2)*e1e2*cross(ev2,jv2/j2t)
     &     -14.*e1j2*j1e2*j1j2*cross(ev2,jv2/j2t)
     &     -7.*e1e2*(8./5.*e1**2-1./5.-7.*e1j2**2
     &     +j1j2**2)*cross(ev2,jv2/j2t))
     &     *oct/TLK*L1/L2/j2t

      j2dot=(j1j2*cross(jv2/j2t,jv)
     &     -5.*e1j2*cross(jv2/j2t,ev))/TLK*L1/L2
     &     -(2.*(e1j2*j1e2*cross(jv2/j2t,jv)
     &     +e1e2*j1j2*cross(jv2/j2t,jv)
     &     +e1j2*j1j2*cross(ev2/e2,jv))
     &     +2.*j1e2*j1j2*cross(jv2/j2t,ev)
     &     -14.*e1e2*e1j2*cross(jv2/j2t,ev)
     &     +(8./5.*e1**2-1./5.-7.*e1j2**2+j1j2**2)*cross(ev2/e2,ev))
     &     *oct/TLK*L1/L2



      yp(7)=e2dot(1)
      yp(8)=e2dot(2)
      yp(9)=e2dot(3)
      yp(10)=j2dot(1)
      yp(11)=j2dot(2)
      yp(12)=j2dot(3)


!     EVOLUTION OF SPINS; SS and SO terms are included
!     vectors S0 and Seff           
       do k=1,3
         S0(k)=(1.+m2/m1)*y(12+k)+(1.+m1/m2)*y(15+k)   
         Seff(k)=(1.+(3./4.)*m2/m1)*y(12+k)
     &        +(1.+(3./4.)*m1/m2)*y(15+k) 
         S1(k)=y(12+k)
         S2(k)=y(15+k)
      end do
      
      S1_SO=2.*L1/(c**2*a1**3*j1**3)*(1.+3./4.*m2/m1)*SP
      S1_SS=m1*m2/(2.*c**2*a1**3*j1**3*(m1+m2)**2)*(1.+m2/m1)*SS
      S2_SO=2.*L1/(c**2*a1**3*j1**3)*(1.+3./4.*m1/m2)*SP
      S2_SS=m1*m2/(2.*c**2*a1**3*j1**3*(m1+m2)**2)*(1.+m1/m2)*SS
      
      S1dot=S1_SO*cross(jv,S1)
     &     +S1_SS*(cross(S0,S1)
     &     -3.*dotp(S0,jv)*cross(jv/j1,S1))    
  
      S2dot=S2_SO*cross(jv,S2)
     &     +S2_SS*(cross(S0,S2)
     &     -3.*dotp(S0,jv)*cross(jv/j1,S2))
      
      yp(13)=S1dot(1)
      yp(14)=S1dot(2)
      yp(15)=S1dot(3)      
      yp(16)=S2dot(1)
      yp(17)=S2dot(2)
      yp(18)=S2dot(3)
            
!     ORBIT reaction
      if(OR.lt.0.5)goto 900
      J1SO=2./(c**2*a1**3*j1**3)*SP
      J1SS=3.*m1*m2/(2.*c**2*L1*a1**3*j1**4*(m1+m2)**2)*SS

      edot=J1SO*(cross(Seff,ev) 
     &     -3.*dotp(Seff,jv/j1)*cross(jv/j1,ev))
     &     +J1SS/2.*(5.*dotp(S0,jv/j1)**2*cross(jv/j1,ev)
     &     -2.*dotp(S0,jv/j1)*cross(S0,ev)
     &     -dotp(S0,S0)*cross(jv/j1,ev))
      jdot=J1SO*cross(Seff,jv)-
     &     J1SS*dotp(S0,jv/j1)*cross(S0,jv)      

      yp(1)=yp(1)+edot(1)      
      yp(2)=yp(2)+edot(2)      
      yp(3)=yp(3)+edot(3)     
      yp(4)=yp(4)+jdot(1)      
      yp(5)=yp(5)+jdot(2)            
      yp(6)=yp(6)+jdot(3)     
 900  continue
      
      return
      end
      
      
C---------------------------------------------------------------------
      SUBROUTINE RK78(IR,T,DT,X,N,TOL,DER)
C---------------------------------------------------------------------
C     Variable step-size automatic one-step integrator for a system of
C     N firts order ordinary differential equations with initial values
C     The Runge-Kutta-Fehlberg formula 7(8) is used
C     REF   E. FEHLBERG, NASA TECHNICAL REPORT TR R-287, 1968
C     Description of parameters list
C     (All floating variables in DOUBLE PRECISION)
C     IR      O    NUMBER OF REJECTIONS OF THE LAST STEP
C     (IN CASE  DT WAS TOO LARGE)
C     T       I-O  INDEPENDENT VARIABLE
C     DT      I-O  STEP SIZE
C     A RECOMMENDED VALUE FOR THE NEXT STEP IS OUTPUT
C     X(N)    I-O  DEPENDENT VARIABLES
C     Fm(N)        AUXILIARY ARRAYS   WITH m = 0 TO 6
C     F7(N)   O    ABSOLUTE ESTIMATED TRUNCATION ERROR ON EACH COMPONENT
C     N       I    ORDER OF THE DIFFERENTIAL EQUATIONS SYSTEM
C     TOL(N)  I    RELATIVE TOLERATED ERROR ON EACH COMPONENT
C     DER     I    NAME OF THE SUBROUTINE COMPUTING THE DERIVATIVES. THIS
C     SUBROUTINE HAS TO HAVE THE STANDARD CALLING SEQUENCE
C     CALL DER(T,X,F0)
C---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      PARAMETER ( CH1 = 34D0/105D0, CH2 = 9D0/35D0, CH3 = 9D0/280D0,
     &     CH4 = 41D0/840D0, AL2 = 2D0/27D0, AL3 = 1D0/9D0,
     &     AL4 = 1D0/6D0,    AL5 = 5D0/12D0, AL6 = 5D-1,
     &     AL7 = 5D0/6D0,    AL9 = 2D0/3D0,  ALA = 1D0/3D0,
     &     B21 = 2D0/27D0,     B31 = 1D0/36D0,    B41 = 1D0/24D0,
     &     B51 = 5D0/12D0,     B61 = 5D-2,        B71 = -25D0/108D0,
     &     B81 = 31D0/3D2,     B101= -91D0/108D0, B111= 2383D0/41D2,
     &     B121= 3D0/205D0,    B131= -1777D0/41D2,B32 = 1D0/12D0,
     &     B43 = .125D0,       B53 = -25D0/16D0,  B64 = 25D-2,
     &     B74 = 125D0/108D0,  B94 = -53D0/6D0,   B104= 23D0/108D0,
     &     B114= -341D0/164D0, B65 = 2D-1,        B75 = -65D0/27D0,
     &     B85 = 61D0/225D0,   B95 = 704D0/45D0,  B105= -976D0/135D0,
     &     B115= 4496D0/1025D0,B76 = 125D0/54D0,  B86 = -2D0/9D0,
     &     B96 = -107D0/9D0,   B106= 311D0/54D0,  B116= -301D0/82D0,
     &     B126= -6D0/41D0,    B136= -289D0/82D0, B87 = 13D0/9D2,
     &     B97 = 67D0/9D1,     B107= -19D0/6D1,   B117= 2133D0/41D2,
     &     B127= -3D0/205D0,   B137= 2193D0/41D2, B108= 17D0/6D0,
     &     B118= 45D0/82D0,    B128= -3D0/41D0,   B138= 51D0/82D0,
     &     B119= 45D0/164D0,   B139= 33D0/164D0,  B1110= 18D0/41D0,
     &     B1310= 12D0/41D0)
C     
      integer IR,n,i
      DIMENSION X(N), TOL(N)
      DIMENSION F0(42),F1(42),F2(42),F3(42),F4(42),F5(42),F6(42),F7(42)
C     

      IF (N .GT. 42) STOP 'N > 42'
C     
      IR = 0
      CALL DER(T, X, F1)
C     
 104  DO I = 1, N
         F0(I) = X(I) + DT*B21*F1(I)
      ENDDO
      CALL DER(T + AL2*DT, F0, F2)
      DO I = 1, N
         F0(I) = X(I) + DT*(B31*F1(I) + B32*F2(I))
      ENDDO
      CALL DER(T + AL3*DT, F0, F3)
      DO I = 1, N
         F0(I) = X(I) + DT*(B41*F1(I) + B43*F3(I))
      ENDDO
      CALL DER(T + AL4*DT, F0, F4)
      DO I = 1, N
         F0(I) = X(I) + DT*(B51*F1(I) + B53*(F3(I) - F4(I)))
      ENDDO
      CALL DER(T + AL5*DT, F0, F5)
      DO I = 1, N
         F0(I) = X(I) + DT*(B61*F1(I) + B64*F4(I) + B65*F5(I))
      ENDDO
      CALL DER(T + AL6*DT, F0, F6)
      DO I = 1, N
         F0(I) = X(I) + DT*(B71*F1(I) + B74*F4(I) + B75*F5(I) +
     &        B76*F6(I))
      ENDDO
      CALL DER(T + AL7*DT, F0, F7)
      DO I = 1, N
         F0(I) = X(I) + DT*(B81*F1(I) + B85*F5(I) + B86*F6(I) +
     &        B87*F7(I))
      ENDDO
      CALL DER(T + AL4*DT, F0, F2)
      DO I = 1, N
         F0(I) = X(I) + DT*(2D0*F1(I) + B94*F4(I) + B95*F5(I) +
     &        B96*F6(I) + B97*F7(I) + 3D0*F2(I))
      ENDDO
      CALL DER(T + AL9*DT, F0, F3)
      DO I = 1, N
         X1 = F1(I)
         X4 = F4(I)
         X5 = F5(I)
         X6 = F6(I)
         X7 = F7(I)
         X8 = F2(I)
         X9 = F3(I)
         F2(I) = CH1*X6 + CH2*(X7 + X8) + CH3*X9
         F0(I) = X(I) + DT*(B101*X1 + B104*X4 + B105*X5 + B106*X6 +
     &        B107*X7 + B108*X8 - B32*X9)
         F4(I) = B111*X1 + B114*X4 + B115*X5 + B116*X6 + B117*X7 +
     &        B118*X8 + B119*X9
         F5(I) = B121*X1 + B126*X6 + B127*X7 + B128*(X8 - X9)
         F6(I) = B131*X1 + B114*X4 + B115*X5 + B136*X6 + B137*X7 +
     &        B138*X8 + B139*X9
      ENDDO
      CALL DER(T + ALA*DT, F0, F3)
      DO I = 1, N
         F7(I) = X(I) + DT*(F4(I) + B1110*F3(I))
         F0(I) = X(I) + DT*(F5(I) - B126*F3(I))
      ENDDO
      CALL DER(T + DT, F7, F4)
      CALL DER(T,      F0, F5)
      DO I = 1, N
         F0(I) = X(I) + DT*(F6(I) + B1310*F3(I) + F5(I))
      ENDDO
      CALL DER(T + DT, F0, F6)
      X7 = 1D-30
      DO I = 1, N
         F0(I) = X(I)
         X(I) = X(I) + DT*(CH3*F3(I) + CH4*(F5(I) + F6(I)) + F2(I))
         F7(I) = DT*(F1(I) + F4(I) - F5(I) - F6(I))*CH4
         X7 = X7 + (F7(I)/TOL(I))**2
      ENDDO
      X9 = DT
      DT = DT*(25D-4/X7)**625D-4
      IF (X7 .GT. 1D0) THEN
         DO I = 1, N
            X(I) = F0(I)
         ENDDO
         IR = IR + 1
         GOTO 104
      ENDIF
      T = T + X9
      RETURN
      END


      end module my_subs





      
      

      
!     MAIN CODE
      
!     Kozai with PN spin evolution
      program kozaispin
      use my_subs
      implicit double precision (a-z)
      integer IR,nq,l,i,j,merge,im,im2,k,h,max,pr
      parameter (nq=19)
      real*8 y(nq),tol(nq),S0(3),Seff(3),q(3),Z1(3),Z2(3)   
      real*8 jv(3),jv2(3),ev(3),ev2(3),jt(3),et(3)
      data tol/nq*1d-16/
      common/cnst/c,pi
      common/orbit/m1,m2,m3,a2,a1
      common/par/GR,SP,SS,OC,OR
      open (unit=68,file="out.dat")

      write(68,*)'#','time,e1(1-3),j1(1-3),e2(1-3),j2(1-3)',
     &     ',S1(1-3),S2(1-3),a1'

      
      pi=2.d0*ASIN(1.d0)
      c=1.d4                    !units M_Sun, AU, G=1
      ut=58./365.               !to convert time in years
      y=0.

      
!     orbital initial conditions
      
      open (unit=99,file="input")
      read(99,*)
      read(99,*)
      read(99,*)e1,e2,a1,a2,M1,M2,M3,chi1,chi2,i1,max,pr
      i1=i1*pi/180.
      
      GR=1.                     !set =0 if you do not want 1pN peri precession, and 2.5pN dissipation
      SP=1.                     !set =0 if you do not want 1pN Spin-Orbit terms
      SS=1.                     !set =0 if you do not want Spin-Spin terms
      OC=1.                     !set =0 if you do not want Octupole terms
      OR=1.                     !set =0 if you do not want nodal LT precession due to spin-orbit terms

      call random(harv)
      an1=0.                    !ascending node
      call random(harv)
      an2=pi-an1
      call random(harv)
      o1=harv*2.*pi   
      call random(harv)
      o2=harv*2.*pi            !harv*2.*pi              
      i2=0.
                           
      
!     INITIAL SPIN VALUES    
      S1=chi1*m1**2/c
      S2=chi2*m2**2/c
      
!     initial vectors
      y(1)=(cos(an1)*cos(o1)-sin(an1)*sin(o1)*cos(i1))*e1 !e1x
      y(2)=(sin(an1)*cos(o1)+cos(an1)*sin(o1)*cos(i1))*e1 !e1y
      y(3)=(sin(o1)*sin(i1))*e1 !e1z

      j1=sqrt(1.-e1**2)
      y(4)=(sin(an1)*sin(i1))*j1 !j1x
      y(5)=(-cos(an1)*sin(i1))*j1 !j1y
      y(6)=(cos(i1))*j1         !j1z
      
      y(7)=(cos(an2)*cos(o2)-sin(an2)*sin(o2)*cos(i2))*e2 !e2x
      y(8)=(sin(an2)*cos(o2)+cos(an2)*sin(o2)*cos(i2))*e2 !e2y
      y(9)=(sin(o2)*sin(i2))*e2 !e2z      

      j2=sqrt(1.-e2**2)
      y(10)=(sin(an2)*sin(i2))*j2 !j2x
      y(11)=(-cos(an2)*sin(i2))*j2 !j2y
      y(12)=(cos(i2))*j2        !j2z      

      
!     initial misalaignment       
      is=i1+0.*pi/180.          !spin-orbit angle

      y(13)=(sin(an1)*sin(is))*S1 !S1x
      y(14)=(-cos(an1)*sin(is))*S1 !S1y
      y(15)=(cos(is))*S1        !S1z

      y(16)=(sin(an1)*sin(is))*S2 !S2x
      y(17)=(-cos(an1)*sin(is))*S2 !S2y
      y(18)=(cos(is))*S2        !S2z

      y(19)=a1
!     initial vectors

      
      Prd=2.*pi*sqrt(a1**3/(m1+m2))
      Pout=2.*pi*sqrt(a2**3/mbh) 
      tau_LK=2.*Prd*(m1+m2)*a2**3*(1.-e2**2)**1.5/(3.*pi*a1**3*m3) !initial Kozai timescale
      DT0=1.e-5                
      DT=DT0
      T=0.       

      i=0
      l=0
      im=0      

      do while(t.lt.max*tau_LK)

!     evolove orbit
         call RK78(IR,t,dt,y,19,tol,force)  
!     evolove orbit
         
         do k=1,3
            S0(k)=(1.+m2/m1)*y(12+k)+(1.+m1/m2)*y(15+k)   
            Seff(k)=(1.+(3./4.)*m2/m1)*y(12+k)
     &           +(1.+(3./4.)*m1/m2)*y(15+k)
            S0(k)=y(12+k)
            Seff(k)=y(15+k)
         end do
         SE=sqrt(Seff(1)**2+Seff(2)**2+Seff(3)**2)
         SO=sqrt(S0(1)**2+S0(2)**2+S0(3)**2)        
         
         do h=1,3
            ev(h)=y(h)
         end do
         do h=1,3
            jv(h)=y(h+3)
         end do
         do h=1,3
            ev2(h)=y(h+6)
         end do
         do h=1,3
            jv2(h)=y(h+9)
         end do
     
         e2=sqrt(dotp(ev2,ev2))
         e1=sqrt(dotp(ev,ev))
         j1=sqrt(dotp(jv,jv))
         j2=sqrt(dotp(jv2,jv2))
         cos_i=dotp(jv,jv2)/abs(j1*j2) !j1*j2
         cost=dotp(Seff,jv2)/abs(SE*j2) !Seff*J2
         cost1=dotp(Seff,jv)/abs(SE*j1) !Seff*J1  
         cosSS=(y(13)*y(16)+y(14)*y(17)+y(15)*y(18))/abs(S1*S2) !S2*S1  
      
         jt=jv*L1+jv2*L2
         jt=jt/sqrt(dotp(jt,jt))

         im=im+1
         if(im.eq.pr)then
            write(68,*)t,y !e1,a1,cos_i,cost,cost1,cosSS
            im=0
         end if
            
         if(y(1).ne.y(1))exit

      end do
      
      end program kozaispin
