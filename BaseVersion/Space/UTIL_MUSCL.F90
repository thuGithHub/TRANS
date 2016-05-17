
REAL FUNCTION ALIMITER(var1,var2)
!!!!!!
	  implicit none
      REAL:: one, small1,var1,var2,epsm
	  
	  one = 1.0
	  small1 = 1.0e-7

        epsm=small1
     
     alimiter=(var1*(var2*var2+2.0*epsm*epsm)+                     &
     &                     var2*(2.0*var1*var1+epsm*epsm))/        &
     &   (2.0*var1*var1-var1*var2+2.0*var2*var2+3.0*epsm*epsm)     &
     &    *0.5*abs(sign(one,var1)+sign(one,var2)) !*0.5 !xiao
	  
	  
END FUNCTION ALIMITER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL FUNCTION MUSCL2(var1,var2)     ! add by dzw,20130817
!!!!!!
	implicit none
    real::var1,var2
      REAL:: one,smooth,epsm
	  
	one = 1.0
    epsm =1.e-6      
   ! var1=2*var1
    !muscl2=sign(one,var2)*max(0.,min(var2*sign(one,var1),2.*var1*sign(one,var2))) 
   smooth=(2.*var1*var2 +epsm)/(var1*var1+var2*var2+epsm)
   muscl2= ((1.+smooth)*var1+ (1.-smooth)*var2)*smooth/2.	  

	  
END FUNCTION MUSCL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


