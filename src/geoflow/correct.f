C*******************************************************************
C     * Copyright (C) 2003 University at Buffalo
C     *
C     * This software can be redistributed free of charge.  See COPYING
C     * file in the top distribution directory for more details.
C     *
C     * This software is distributed in the hope that it will be useful,
C     * but WITHOUT ANY WARRANTY; without even the implied warranty of
C     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C     *
C     * Author: 
C     * Description: 
C     *
C*******************************************************************
C     * $Id: correct.f 143 2007-06-25 17:58:08Z dkumar $ 
C     *

C***********************************************************************
      subroutine correct(uvec, uprev, fluxxp, fluxyp, fluxxm, fluxym,
     1     tiny, dtdx, dtdy, dt, dUdx, dUdy,lap_phi,
     2     curv,intfrictang, bedfrictang, g, kactxy,  dgdx, 
     3     frict_tiny, forceint,forcebed, dragfoce ,DO_EROSION,
     4     eroded, Vel, terminal_vel,
     5     eps, IF_STOPPED, fluxsrc,eta,min_dx,elsrelti)
C***********************************************************************

      implicit none
      double precision forceint, forcebed, eroded, speed
      double precision forceintx, forceinty
      double precision forcebedx, forcebedy
      double precision forcebedmax, forcebedequil, forcegrav
      double precision unitvx, unitvy, Vel(2)
      double precision alphaxx, alphayy, alphaxy, alphayz
      double precision tanbed, terminal_vel, dragfoce(2)
      double precision fluxxp(6),fluxyp(6),tiny, uprev(6), ustore(6)
      double precision fluxxm(6), fluxym(6)
      double precision uvec(6), dUdx(6), dUdy(6),lap_phi(2)
      double precision h_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx,tmp
      double precision dnorm, fluxsrc(6)
      double precision t1, t2, t3, t4,dvx,dvy
      double precision erosion_rate,threshold,es,gama
      double precision eps, drag(4),eta,elsrelti,min_dx,cap_width

!     function calls
      double precision sgn

      integer i
      integer DO_EROSION, IF_STOPPED
      parameter(threshold=1.0D-02,erosion_rate=0.1)

c     initialize to zero
      forceintx=0.d0
      forcebedx=0.d0
      forceinty=0.d0
      forcebedy=0.d0
      unitvx=0.d0
      unitvy=0.d0
      eroded=0.d0

c     -------------------------------Hossein-------------------------------------

      do 10 i = 2,4
         ustore(i)=uprev(i)+dt*fluxsrc(i)
     1        -dtdx*(fluxxp(i)-fluxxm(i))
     2        -dtdy*(fluxyp(i)-fluxym(i))
 10   continue

      ustore(1)=uprev(1)+dt*fluxsrc(1)
     $     -dtdx*(fluxxp(1)+fluxxm(5))
     $     -dtdy*(fluxyp(1)+fluxym(5))

      cap_width=6*min_dx
      gama=.01*min_dx*min_dx

c     we set relaxation_time to (min_dx**2)
c     and eta to 5 min_dx so 1/eta**2=1/(25*(min_dx**2))
c     2 terms of RHS of phase field equation is updated in explicit
c     time scheme which are t2 or lagrange multiplier
c     to conserve the mass and d_F/d_phi term or t1
c     for t1 we have dt*((min_dx**2)/dt*eta*d_F/d_phi)
c     which is equal to t1=.04*uvec(1)*(uvec(1)**2-1)
      t1=gama*(uvec(1)*(uvec(1)**2-1)/(cap_width*cap_width))

      t2=gama*eta*(uvec(1)**2-1)/(cap_width*cap_width)


      ustore(1)=ustore(1)+dt*(t2-t1)
c      (f-eta)/(capwid^2)  capwid=4*mindx
!      print *,"coefficient is:  ", dt*elsrelti/(16*capwid*capwid)
!      print *,"source is:       ", -t3+eta
      
      ustore(2) = max(ustore(2),0.d0)
c      ustore(1) = max(ustore(1),0.)
      ustore(5) = uvec(5)
      ustore(6) = uvec(6)

      if(uvec(2).gt.tiny) then !tiny
c     Source terms ...
c     here speed is speed squared
         speed=Vel(1)**2+Vel(2)**2
         if(speed.gt.0.d0) then
c     here speed is speed
            speed=dsqrt(speed)
            unitvx=Vel(1)/speed
            unitvy=Vel(2)/speed
         else
            unitvx=0.d0
            unitvy=0.d0
         endif
         tanbed=dtan(bedfrictang)
         h_inv = 1.0/uvec(2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     x-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     alphaxy -- see pitman-le (2005)
         dvx = h_inv*(dUdy(3)-Vel(1)*dUdy(2))
         sgn_dudy = sgn(dvx, frict_tiny)
         alphaxy = sgn_dudy*dsin(intfrictang)*kactxy

         t2=alphaxy*uvec(2)*(g(3)*dUdy(2)
     $        +dgdx(2)*uvec(2))
         t3=unitvx*
     $        dmax1(g(3)*uvec(2)+Vel(1)*uvec(3)*curv(1),0.d0)
     $        *tanbed

         t4 = uvec(2)*g(1)


c     update ustore
         ustore(3) = ustore(3) + dt*(t4-t2-t3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     solid fraction y-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         dvy = h_inv*(dUdx(4)-Vel(2)*dUdx(2))
         sgn_dvdx = sgn(dvy, frict_tiny)
         alphaxy = sgn_dvdx*dsin(intfrictang)*kactxy

         t2=alphaxy*uvec(2)*(g(3)*dUdx(2)
     $        +dgdx(1)*uvec(2))
c     ------------------------------------------------  the internal friction force ------------------------------------
         t3=unitvy*
     $        dmax1(g(3)*uvec(2)+Vel(2)*uvec(4)*curv(2),0.d0)
     $        *tanbed

c-------------------------------the bed friction force for fast moving flow---------------------------------------
         t4 = uvec(2)*g(2)

         ustore(4) = ustore(4) + dt*(t4-t2-t3)
      endif

      do 20 i = 1,4
         uvec(i)=ustore(i)
 20   continue

      return
      end
