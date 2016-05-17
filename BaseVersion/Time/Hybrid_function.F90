	 REAL function FSST(Den,vK,vOmega,Dist,RmiuL,CDiff)
!!!!!!
    USE Global
    IMPLICIT NONE
!!!!!!
     REAL::Den,vK,vOmega,Dist,RmiuL,CDiff
     REAL::arg1,arg2,arga,argb,temp,arg
!	  Blending Function
!!!!!!
	  if(Kind_subHY == 0)then													!!original Strelts' DES!!
	     Fsst=0.

	  elseif(Kind_subHY == 1)then												!! Menter-DDES-F1
		arg1=sqrt(Vk)/(0.09*Vomega*DIST*Ref)
		arg2=500.*RmiuL/(Den*DIST*DIST*Vomega*Ref*Ref)
		arga=max(arg1,arg2)
		temp=max(Den*CDiff*Ref,1.D-10)										!! CD|k-w cross derivative term 
		argb=4.*Den*0.856*Vk/(temp*DIST*DIST)
		arg =min(arga,argb)
		Fsst=tanh(arg*arg*arg*arg)

	  elseif (Kind_subHY == 2 )	then											!! Menter-DDES-F2
		arg1=2.*sqrt(Vk)/(0.09*Vomega*DIST*Ref)
		arg2=500.*RmiuL/(Den*DIST*DIST*Vomega*Ref*Ref)
		arg=max(arg1,arg2)
		Fsst=tanh(arg*arg)

	  elseif (Kind_subHY == 3 ) then											!! Spalart-DDES-F1
		arg1=sqrt(Vk)/(0.09*Vomega*DIST*Ref)
		arg2=500.*RmiuL/(Den*DIST*DIST*Vomega*Ref*Ref)
		arga=max(arg1,arg2)
		temp=max(Den*CDiff*Ref,1.D-10)										!! CD|k-w cross derivative term 
		argb=4.*Den*0.856*Vk/(temp*DIST*DIST)
		arg =min(arga,argb)
		Fsst=tanh(arg*arg*arg*arg)
	  endif
!
	  return
	  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  !!!求解IDDES中的两个函数
	 REAL function  fd_function (Den1,   RmiuT1,  aDst1,  aLmax1,dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz             )
	   
    USE Global
    IMPLICIT NONE
	 	REAL:: Den1,   RmiuT1,  aDst1,  aLmax1,dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
		REAL:: alpha1,Sij_square,rdt,fdt,fb,fd_real  
		 
		 alpha1=0.25-aDst1/aLmax1
		 fb=min( 2.* exp(-9.0*alpha1*alpha1) , 1. )
           
           Sij_square = dudx*dudx + dudy*dudy +  dudz*dudz + dvdx*dvdx + dvdy*dvdy +  dvdz*dvdz + dwdx*dwdx + dwdy*dwdy +  dwdz*dwdz
      	 
		 Sij_square = sqrt(Sij_square)
		 Sij_square = max(Sij_square,1.D-10)  
		         
		 rdt=(RmiuT1/Den1)/(akapa_iddes**2.* aDst1**2.*Sij_square*Ref)
		 fdt=1.-tanh((8.* rdt)*(8.* rdt)*(8.* rdt))
		 
		 fd_function= max( (1.-fdt), fb)
         fd_real = fb/(1.-fdt)
	    
	  return
	  end
	  
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  
    REAL function  fe_function (Den2,   RmiuT2,  RMIUL2,   aDst2,  aLmax2, dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz      ) 
	 
    USE Global
    IMPLICIT NONE
	  	REAL::Den2,   RmiuT2,  RMIUL2,   aDst2,  aLmax2, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
		REAL::alpha2,fe1, Sij_square,fe2,rdt,ft,rdl,fl
		
		alpha2=0.25-aDst2/aLmax2
          
		if ( alpha2 >= 0. ) then		 
		 fe1=2.* exp(-11.09 * alpha2*alpha2)	    
		else
		 fe1=2.* exp(-9. * alpha2*alpha2)
		end if
		
           Sij_square = dudx*dudx + dudy*dudy +  dudz*dudz+ dvdx*dvdx + dvdy*dvdy +  dvdz*dvdz+ dwdx*dwdx + dwdy*dwdy +  dwdz*dwdz
      	 
		Sij_square = sqrt(Sij_square)
		Sij_square = max(Sij_square,1.D-10)  
		         
		rdt=(RmiuT2/Den2)/(akapa_iddes**2.* aDst2**2.*Sij_square*Ref)
	    ft=tanh( (Ct**2.*rdt) * (Ct**2.*rdt) * (Ct**2.*rdt) )



		         
		rdl=(RMIUL2/Den2)/(akapa_iddes**2.* aDst2**2.*Sij_square*Ref)
        fl= tanh( (Cl**2.*rdl)**10. )
		
		fe2 = 1. - max( ft, fl)

		fe_function = fe2 * max( (fe1-1.), 0. )

	  
	  return
	  end  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!