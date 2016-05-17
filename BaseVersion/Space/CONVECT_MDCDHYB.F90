SUBROUTINE CONVECT_MDCDHYB   !by dzw05,20130529

    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,P,M,N
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER::  npt

     REAL:: LU(6,5),RU(6,5) 
     REAL:: ISK(4),ak(4)
	 REAL:: v_3, v_2, v_1, v0, v1, v2, v3  ! xprime
	 REAL:: w0,w1,w2,w3
     REAL:: aw0,aw1,aw2,aw3                 !gama_min

	 REAL:: qq0, qq1, qq2,qq3,beta
     REAL:: epsm,sensor,gama

   INTERFACE
        REAL FUNCTION alimiter(var1,var2)
        END FUNCTION
    END INTERFACE
      
       beta=sqrt(13.0/12.0)
       epsm=1.0d-6   


    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1







  !    I direction
  !--------------------------------------------------------------    
 do K=1,NK1
 do J=1,NJ1
 do i=1,NI        
      do npt=1,6   
       LU(npt,1:5)=(/V(2,i-4+npt,j,k)/V(1,i-4+npt,j,k),V(3,i-4+npt,j,k)/V(1,i-4+npt,j,k),V(4,i-4+npt,j,k)/V(1,i-4+npt,j,k),PP(i-4+npt,j,k),T(i-4+npt,j,k)/) 
      enddo
     
     do npt=1,6       
       RU(npt,1:5)=(/V(2,i-4+npt,j,k)/V(1,i-4+npt,j,k),V(3,i-4+npt,j,k)/V(1,i-4+npt,j,k),V(4,i-4+npt,j,k)/V(1,i-4+npt,j,k),PP(i-4+npt,j,k),T(i-4+npt,j,k)/)
      enddo

 !!==================================================================
 !! Ql(i)
 !!==================================================================
 
    do L=1,5
	  v_2= LU(1,L)
	  v_1= LU(2,L)
      v0 = LU(3,L)
	  v1 = LU(4,L)
	  v2 = LU(5,L)
	  v3 = LU(6,L) 
	 
	  
      sensor=max(shock(i-3,j,k),shock(i-2,j,k),shock(i-1,j,k),shock(i,j,k),shock(i+1,j,k),shock(i+2,j,k))
      !sensor =0.0
      !do P=i-3,i+2
	  ! do M=j-1,j+1
	  !  do N=k-1,k+1
     !	sensor=max(sensor,Shock(P,M,N))
     !   enddo
	 ! enddo
	!enddo
            
       gama=gama_min+sensor*(gama_max-gama_min)   
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)

      !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
      
      ISK(1)=( 0.5*v_2 -2.0*v_1 +1.5*v0)**2+(beta*v_2  -2.0*beta*v_1 +beta*v0)**2
	  ISK(2)=(-0.5*v_1 +0.0*v0  +0.5*v1)**2+(beta*v_1  -2.0*beta*v0  +beta*v1)**2
	  ISK(3)=(-1.5*v0  +2.0*v1  -0.5*v2)**2+(beta*v0   -2.0*beta*v1  +beta*v2)**2
      ISK(4)= (2.5*v1  -4.0*v2  +1.5*v3)**2+(beta*v1   -2.0*beta*v2  +beta*v3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif

      qq0=  2./6.*v_2 - 7./6.*v_1 + 11./6.*v0
	  qq1=- 1./6.*v_1 + 5./6.*v0  +  2./6.*v1
	  qq2=  2./6.*v0  + 5./6.*v1  -  1./6.*v2
      qq3= 11./6.*v1   -7./6.*v2  +  1./3.*v3
	
    VL(1,L,i,j,k)=w0*qq0+w1*qq1+w2*qq2+w3*qq3

    if(sensor>= F_limiter2 )then
         VL(1,L,i,j,k)=v0 + 0.5*alimiter(v1-v0,v0-v_1)
    endif
 
!!===============================
 !! Qr(i)
 !!==============================
             
     v_3= RU(1,L)
     v_2= RU(2,L)
	 v_1= RU(3,L)
     v0 = RU(4,L)
	  v1= RU(5,L)
	  v2= RU(6,L)
    
!      sensor=max(shock(i-3,j,k),shock(i-2,j,k),shock(i-1,j,k),shock(i,j,k),shock(i+1,j,k),shock(i+2,j,k))
!       gama=gama_min+sensor*(0.0463783-gama_min)   
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)
      
      
          !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
      
      ISK(1)=( 0.5*v2 -2.0*v1 +1.5*v0)**2+(beta*v2  -2.0*beta*v1 +beta*v0)**2
	  ISK(2)=(-0.5*v1 +0.0*v0  +0.5*v_1)**2+(beta*v1  -2.0*beta*v0  +beta*v_1)**2
	  ISK(3)=(-1.5*v0  +2.0*v_1  -0.5*v_2)**2+(beta*v0   -2.0*beta*v_1  +beta*v_2)**2
      ISK(4)= (2.5*v_1  -4.0*v_2  +1.5*v_3)**2+(beta*v_1   -2.0*beta*v_2  +beta*v_3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif


      qq0=  2./6.*v2  - 7./6.*v1  + 11./6.*v0  
	  qq1= -1./6.*v1  + 5./6.*v0  +  2./6.*v_1
	  qq2=  2./6.*v0  + 5./6.*v_1 -  1./6.*v_2
      qq3= 11./6.*v_1 - 7./6.*v_2  + 1./3.*v_3

     VR(1,L,i,j,k)=w0*qq0 +w1*qq1  +w2*qq2 +w3*qq3    
    
    if(sensor>= F_limiter2 )then
         VR(1,L,i,j,k)=v0 - 0.5*alimiter(v0-v_1,v1-v0)
    endif	  
   enddo
enddo
enddo
enddo



 !    J direction
  !--------------------------------------------------------------    
 do K=1,NK1
 do i=1,NI1  
 do J=1,NJ
      
      do npt=1,6   
       LU(npt,1:5)=(/V(2,i,j-4+npt,k)/V(1,i,j-4+npt,k),V(3,i,j-4+npt,k)/V(1,i,j-4+npt,k),V(4,i,j-4+npt,k)/V(1,i,j-4+npt,k),PP(i,j-4+npt,k),T(i,j-4+npt,k)/) 
      enddo
     
     do npt=1,6       
       RU(npt,1:5)=(/V(2,i,j-4+npt,k)/V(1,i,j-4+npt,k),V(3,i,j-4+npt,k)/V(1,i,j-4+npt,k),V(4,i,j-4+npt,k)/V(1,i,j-4+npt,k),PP(i,j-4+npt,k),T(i,j-4+npt,k)/) 
      enddo

 !!==================================================================
 !! Ql(j)
 !!==================================================================
 
    do L=1,5
	  v_2= LU(1,L)
	  v_1= LU(2,L)
      v0 = LU(3,L)
	  v1 = LU(4,L)
	  v2 = LU(5,L)
	  v3 = LU(6,L) 
	 
	  sensor=max(shock(i,j-3,k),shock(i,j-2,k),shock(i,j-1,k),shock(i,j,k),shock(i,j+1,k),shock(i,j+2,k))
!       sensor =0.0
!      do P=i-1,i+1
!	   do M=j-3,j+2
!	    do N=k-1,k+1
 !    	sensor=max(sensor,Shock(P,M,N))
  !      enddo
	!  enddo
	!enddo
    
        gama=gama_min+sensor*(gama_max-gama_min)   
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)
            !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
      
      ISK(1)=( 0.5*v_2 -2.0*v_1 +1.5*v0)**2+(beta*v_2  -2.0*beta*v_1 +beta*v0)**2
	  ISK(2)=(-0.5*v_1 +0.0*v0  +0.5*v1)**2+(beta*v_1  -2.0*beta*v0  +beta*v1)**2
	  ISK(3)=(-1.5*v0  +2.0*v1  -0.5*v2)**2+(beta*v0   -2.0*beta*v1  +beta*v2)**2
      ISK(4)= (2.5*v1  -4.0*v2  +1.5*v3)**2+(beta*v1   -2.0*beta*v2  +beta*v3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif
      

      qq0=  2./6.*v_2 - 7./6.*v_1 + 11./6.*v0
	  qq1=- 1./6.*v_1 + 5./6.*v0  +  2./6.*v1
	  qq2=  2./6.*v0  + 5./6.*v1  -  1./6.*v2
      qq3= 11./6.*v1   -7./6.*v2  +  1./3.*v3
	
    VL(2,L,i,j,k)=w0*qq0+w1*qq1+w2*qq2+w3*qq3

    if(sensor>= F_limiter2 )then
         VL(2,L,i,j,k)=v0 + 0.5*alimiter(v1-v0,v0-v_1)
    endif
 
!!===============================
 !! Qr(j)
 !!==============================
             
     v_3=RU(1,L)
     v_2=RU(2,L)
	 v_1=RU(3,L)
     v0=RU(4,L)
	  v1= RU(5,L)
	  v2= RU(6,L)
    
 !    sensor=max(shock(i,j-3,k),shock(i,j-2,k),shock(i,j-1,k),shock(i,j,k),shock(i,j+1,k),shock(i,j+2,k))
  !     gama=gama_min+sensor*(0.0463783-gama_min)   
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)
      
           !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
      
      ISK(1)=( 0.5*v2 -2.0*v1 +1.5*v0)**2+(beta*v2  -2.0*beta*v1 +beta*v0)**2
	  ISK(2)=(-0.5*v1 +0.0*v0  +0.5*v_1)**2+(beta*v1  -2.0*beta*v0  +beta*v_1)**2
	  ISK(3)=(-1.5*v0  +2.0*v_1  -0.5*v_2)**2+(beta*v0   -2.0*beta*v_1  +beta*v_2)**2
      ISK(4)= (2.5*v_1  -4.0*v_2  +1.5*v_3)**2+(beta*v_1   -2.0*beta*v_2  +beta*v_3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif


      qq0=  2./6.*v2  - 7./6.*v1  + 11./6.*v0  
	  qq1= -1./6.*v1  + 5./6.*v0  +  2./6.*v_1
	  qq2=  2./6.*v0  + 5./6.*v_1 -  1./6.*v_2
      qq3= 11./6.*v_1 - 7./6.*v_2  + 1./3.*v_3

VR(2,L,i,j,k)=w0*qq0 +w1*qq1  +w2*qq2 +w3*qq3 
    if(sensor>= F_limiter2 )then
         VR(2,L,i,j,k)=v0 - 0.5*alimiter(v0-v_1,v1-v0)
    endif	   
	  
   enddo
enddo
enddo
enddo


 !    K direction
  !--------------------------------------------------------------    

 do J=1,NJ1
 do i=1,NI1  
 do K=1,NK      
      do npt=1,6   
       LU(npt,1:5)=(/V(2,i,j,k-4+npt)/V(1,i,j,k-4+npt),V(3,i,j,k-4+npt)/V(1,i,j,k-4+npt),V(4,i,j,k-4+npt)/V(1,i,j,k-4+npt),PP(i,j,k-4+npt),T(i,j,k-4+npt)/) 
      enddo
     
     do npt=1,6       
       RU(npt,1:5)=(/V(2,i,j,k-4+npt)/V(1,i,j,k-4+npt),V(3,i,j,k-4+npt)/V(1,i,j,k-4+npt),V(4,i,j,k-4+npt)/V(1,i,j,k-4+npt),PP(i,j,k-4+npt),T(i,j,k-4+npt)/)  
      enddo

 !!==================================================================
 !! Ql(k)
 !!==================================================================
 
    do L=1,5
	  v_2= LU(1,L)
	  v_1= LU(2,L)
      v0 = LU(3,L)
	  v1 = LU(4,L)
	  v2 = LU(5,L)
	  v3 = LU(6,L) 
	 
	 sensor=max(shock(i,j,k-3),shock(i,j,k-2),shock(i,j,k-1),shock(i,j,k),shock(i,j,k+1),shock(i,j,k+2))
!       sensor =0.0
 !     do P=i-1,i+1
!	   do M=j-1,j+1
!	    do N=k-3,k+2
!     	sensor=max(sensor,Shock(P,M,N))
!        enddo
!	  enddo
!	enddo	
        
       gama=gama_min+sensor*(gama_max-gama_min) 
     !  write(*,*)sensor
      ! pause  
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)
      
          !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
      
      ISK(1)=( 0.5*v_2 -2.0*v_1 +1.5*v0)**2+(beta*v_2  -2.0*beta*v_1 +beta*v0)**2
	  ISK(2)=(-0.5*v_1 +0.0*v0  +0.5*v1)**2+(beta*v_1  -2.0*beta*v0  +beta*v1)**2
	  ISK(3)=(-1.5*v0  +2.0*v1  -0.5*v2)**2+(beta*v0   -2.0*beta*v1  +beta*v2)**2
      ISK(4)= (2.5*v1  -4.0*v2  +1.5*v3)**2+(beta*v1   -2.0*beta*v2  +beta*v3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif

      qq0=  2./6.*v_2 - 7./6.*v_1 + 11./6.*v0
	  qq1=- 1./6.*v_1 + 5./6.*v0  +  2./6.*v1
	  qq2=  2./6.*v0  + 5./6.*v1  -  1./6.*v2
      qq3= 11./6.*v1   -7./6.*v2  +  1./3.*v3
	
    VL(3,L,i,j,k)=w0*qq0+w1*qq1+w2*qq2+w3*qq3

    if(sensor>= F_limiter2 )then
         VL(3,L,i,j,k)=v0 + 0.5*alimiter(v1-v0,v0-v_1)
    endif
 
!!===============================
 !! Qr(k)
 !!==============================
             
     v_3=RU(1,L)
     v_2=RU(2,L)
	 v_1=RU(3,L)
     v0=RU(4,L)
	  v1= RU(5,L)
	  v2= RU(6,L)
    
 !     sensor=max(shock(i,j,k-3),shock(i,j,k-2),shock(i,j,k-1),shock(i,j,k),shock(i,j,k+1),shock(i,j,k+2))
        
        
  !     gama=gama_min+sensor*(0.0463783-gama_min)   
       
       aw0=1.5*(0.0463783+gama)
	   aw1=0.5-1.5*0.0463783+4.5*gama
	   aw2=0.5-1.5*0.0463783-4.5*gama
       aw3=1.5*(0.0463783-gama)
      
          !if (sensor < epsm) then

     ! w0=aw0
     ! w1=aw1
     ! w2=aw2
      !w3=aw3
     ! else
  
      ISK(1)=( 0.5*v2 -2.0*v1 +1.5*v0)**2+(beta*v2  -2.0*beta*v1 +beta*v0)**2
	  ISK(2)=(-0.5*v1 +0.0*v0  +0.5*v_1)**2+(beta*v1  -2.0*beta*v0  +beta*v_1)**2
	  ISK(3)=(-1.5*v0  +2.0*v_1  -0.5*v_2)**2+(beta*v0   -2.0*beta*v_1  +beta*v_2)**2
      ISK(4)= (2.5*v_1  -4.0*v_2  +1.5*v_3)**2+(beta*v_1   -2.0*beta*v_2  +beta*v_3)**2
      ISK(4)=max(ISK(1),ISK(2),ISK(3),ISK(4))
    
      ak(1)=aw0/(epsm+ISK(1))**2
	  ak(2)=aw1/(epsm+ISK(2))**2
	  ak(3)=aw2/(epsm+ISK(3))**2 
      ak(4)=aw3/(epsm+ISK(4))**2 
        
	  w0=ak(1)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w1=ak(2)/(ak(1)+ak(2)+ak(3)+ak(4))
	  w2=ak(3)/(ak(1)+ak(2)+ak(3)+ak(4))
      w3=ak(4)/(ak(1)+ak(2)+ak(3)+ak(4))
     ! endif


      qq0=  2./6.*v2  - 7./6.*v1  + 11./6.*v0  
	  qq1= -1./6.*v1  + 5./6.*v0  +  2./6.*v_1
	  qq2=  2./6.*v0  + 5./6.*v_1 -  1./6.*v_2
      qq3= 11./6.*v_1 - 7./6.*v_2  + 1./3.*v_3

VR(3,L,i,j,k)=w0*qq0 +w1*qq1  +w2*qq2 +w3*qq3 
   
    if(sensor>= F_limiter2 )then
         VR(3,L,i,j,k)=v0 - 0.5*alimiter(v0-v_1,v1-v0)
    endif	   
	  
   enddo
enddo
enddo
enddo









END SUBROUTINE CONVECT_MDCDHYB
