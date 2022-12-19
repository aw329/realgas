c--------------------------------------------------------------------	  
      subroutine getpt(dgas,ugas,tgas,pgas)
c--------------------------------------------------------------------	        
c dgas = density (kg/m^3)
c ugas = internal energy (J/kg)
c tgas = temperature (K)
c pgas = pressure (Pa)
c zgas = compressibility factor
c dref = reference density (kg/m^3)
c rconst = gas constant (J/kg/K)
c rrconst = 1/rconst
c--------------------------------------------------------------------
c use module containing coefficients
      use realgas
c      
      implicit none  
      integer m,n
      double precision 
     & tt,fn,fnm1,fm,logd,dd,logd2,logd3,logd4,
     & d2,d3,ygasnm1,ygasn,dgas,ugas,tgas,pgas,
     & xgas,ygas,zgas,xgas2,xgas3  
c
      logd = dlog(dgas/dref)  
      dd = dgas-dref
      d2 = dgas*dgas - dref*dref
      d3 = dgas*dgas*dgas - dref*dref*dref
c      
      xgas = (dgas-dm)/df
      ygas = (ugas-um)/uf
      xgas2 = xgas*xgas 
      xgas3 = xgas2*xgas 
      logd2 = (dd - dm*logd)*rdf
      logd3 = (dm2*logd + d2 - 4.0*dm*dd)*rdf2
      logd4 =-(dm3*logd - 1.5d0*dm2*dd + 1.5d0*dm*d2 - third*d3)*rdf3   
c
      if(dgas.ge.dref) then
      tt = 0.0d0
      zgas = PHI(1,1) + PHI(2,1)*xgas + PHI(3,1)*xgas2 + PHI(4,1)*xgas3
      zgas = zgas + PHI(1,2)*ygas  
      tt = tt - PHI(1,2)*logd 
      zgas = zgas + PHI(2,2)*xgas*ygas	  
      tt = tt -  PHI(2,2)*logd2
      zgas = zgas + PHI(3,2)*xgas2*ygas	  	  
      tt = tt -  PHI(3,2)*logd3
      zgas = zgas + PHI(4,2)*xgas3*ygas	  	  
      tt = tt -  PHI(4,2)*logd4 + aid(5)	  	 
      ygasnm1 = 2.0d0*ygas	  
      ygasn = ygas*ygas	  
      zgas = zgas + PHI(1,3)*ygasn  
      tt = tt - PHI(1,3)*ygasnm1*logd 
      zgas = zgas + PHI(2,3)*xgas*ygasn	  
      tt = tt -  PHI(2,3)*ygasnm1*logd2
      zgas = zgas + PHI(3,3)*xgas2*ygasn	  	  
      tt = tt -  PHI(3,3)*ygasnm1*logd3
      zgas = zgas + PHI(4,3)*xgas3*ygasn	  	  
      tt = tt -  PHI(4,3)*ygasnm1*logd4 + aid(4)*ygasnm1
      ygasnm1 = 3.0*ygasn	  
      ygasn = ygasn*ygas	  
      zgas = zgas + PHI(1,4)*ygasn  
      tt = tt - PHI(1,4)*ygasnm1*logd 
      zgas = zgas + PHI(2,4)*xgas*ygasn	  
      tt = tt -  PHI(2,4)*ygasnm1*logd2
      zgas = zgas + PHI(3,4)*xgas2*ygasn	  	  
      tt = tt -  PHI(3,4)*ygasnm1*logd3
      zgas = zgas + PHI(4,4)*xgas3*ygasn	  	  
      tt = tt -  PHI(4,4)*ygasnm1*logd4 + aid(3)*ygasnm1	  
      ygasnm1 = 4.0d0*ygasn	  
      ygasn = ygasn*ygas	  
      zgas = zgas + PHI(1,5)*ygasn  
      tt = tt - PHI(1,5)*ygasnm1*logd 
      zgas = zgas + PHI(2,5)*xgas*ygasn	  
      tt = tt -  PHI(2,5)*ygasnm1*logd2
      zgas = zgas + PHI(3,5)*xgas2*ygasn	  	  
      tt = tt -  PHI(3,5)*ygasnm1*logd3
      zgas = zgas + PHI(4,5)*xgas3*ygasn	  	  
      tt = tt -  PHI(4,5)*ygasnm1*logd4 + aid(2)*ygasnm1
      ygasnm1 = 5.0d0*ygasn	  
      ygasn = ygasn*ygas	  
      zgas = zgas + PHI(1,6)*ygasn  
      tt = tt - PHI(1,6)*ygasnm1*logd 
      zgas = zgas + PHI(2,6)*xgas*ygasn	  
      tt = tt -  PHI(2,6)*ygasnm1*logd2
      zgas = zgas + PHI(3,6)*xgas2*ygasn	  	  
      tt = tt -  PHI(3,6)*ygasnm1*logd3
      zgas = zgas + PHI(4,6)*xgas3*ygasn	  	  
      tt = tt -  PHI(4,6)*ygasnm1*logd4 + aid(1)*ygasnm1
      else	 
      tt =       5.0d0*aid(1)*(ygas**4.0d0)
     &         + 4.0d0*aid(2)*(ygas**3.0d0)
     &         + 3.0d0*aid(3)*(ygas**2.0d0) 
     &         + 2.0d0*aid(4)*ygas + aid(5)
      zgas = 1.0d0
      endif 
c
      if(dgas.lt.0.01d0) dgas = 0.01d0	  	  
c  
      tgas = uf/tt*rrconst
      if(tgas.lt.188.0d0) tgas = 188.0d0	  
      pgas = rconst*zgas*tgas*dgas
      if(pgas.lt.0.01d0) then 
      pgas = 0.01d0  
      dgas = pgas/(rconst*tgas)
      endif	
c	  
      return
      end	  
c
c--------------------------------------------------------------------	
