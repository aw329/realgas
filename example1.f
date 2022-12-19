c--------------------------------------------------------------------	  
      subroutine gets(dgas,ugas,ds)
c--------------------------------------------------------------------	  
c ds = entropy (J/kg/K)     
c dgas = density (kg/m^3)
c ugas = internal energy (J/kg)
c dref = reference density (kg/m^3)
c rconst = gas constant (J/kg/K)
c rrconst = 1/rconst
c PHI contains alpha coefficients for compressibility factor 
c PHI(1,1)=alpha_00 ... etc
c aid contains beta coefficients 
c beta_0 = aid(5), beta_1 = aid(4) ... etc
c--------------------------------------------------------------------
c module containing coefficients
      use realgas
c   
      implicit none  
c	  
      integer m,n
      double precision ds,fn,fnm1,fm,logd,dd,logd2,logd3,
     & logd4,d2,d3,ygasn,dgas,ugas,xgas,ygas	  
c	  
      logd = dlog(dgas/dref)	  
      dd = dgas-dref
      d2 = dgas*dgas - dref*dref
      d3 = dgas*dgas*dgas - dref*dref*dref
      xgas = (dgas-dm)/df	
      ygas = (ugas-um)/uf	
      logd2 = (dd - dm*logd)*rdf
      logd3 = (dm2*logd + d2 - 4.0d0*dm*dd)*rdf2
      logd4 =-(dm3*logd - 1.5d0*dm2*dd + 1.5d0*dm*d2 - third*d3)*rdf3  
c
      if(dgas.ge.dref) then
      ds = 0.0d0
      ds = ds - PHI(1,1)*logd 
      ds = ds - PHI(2,1)*logd2	  
      ds = ds - PHI(3,1)*logd3	 
      ds = ds - PHI(4,1)*logd4	+ aid(6)       
      ds = ds - PHI(1,2)*logd*ygas 
      ds = ds - PHI(2,2)*logd2*ygas	  
      ds = ds - PHI(3,2)*logd3*ygas  
      ds = ds - PHI(4,2)*logd4*ygas + aid(5)*ygas 
      ygasn = ygas*ygas	  
      ds = ds - PHI(1,3)*logd*ygasn 
      ds = ds - PHI(2,3)*logd2*ygasn	  
      ds = ds - PHI(3,3)*logd3*ygasn 
      ds = ds - PHI(4,3)*logd4*ygasn + aid(4)*ygasn	       
      ygasn = ygasn*ygas	
      ds = ds - PHI(1,4)*logd*ygasn 
      ds = ds - PHI(2,4)*logd2*ygasn	  
      ds = ds - PHI(3,4)*logd3*ygasn 
      ds = ds - PHI(4,4)*logd4*ygasn + aid(3)*ygasn	       
      ygasn = ygasn*ygas	
      ds = ds - PHI(1,5)*logd*ygasn 
      ds = ds - PHI(2,5)*logd2*ygasn	  
      ds = ds - PHI(3,5)*logd3*ygasn 
      ds = ds - PHI(4,5)*logd4*ygasn + aid(2)*ygasn	  
      ygasn = ygasn*ygas	
      ds = ds - PHI(1,6)*logd*ygasn 
      ds = ds - PHI(2,6)*logd2*ygasn	  
      ds = ds - PHI(3,6)*logd3*ygasn 	  
      ds = ds - PHI(4,6)*logd4*ygasn + aid(1)*ygasn	        
      else	 
      ds = aid(1)*(ygas**5.0d0)
     & +   aid(2)*(ygas**4.0d0)
     & +   aid(3)*(ygas**3.0d0) 
     & +   aid(4)*(ygas*ygas) 
     & +   aid(5)*ygas + aid(6)
      endif 
c
      ds = ds*rconst  
c	  
      return
      end 
c
c--------------------------------------------------------------------	
