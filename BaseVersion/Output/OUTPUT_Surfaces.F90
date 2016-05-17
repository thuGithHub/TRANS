SUBROUTINE OUTPUT_Surfaces
    USE Global
    IMPLICIT NONE
    
    INTEGER:: iBlock
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        CALL OUTPUT_Surface
    ENDDO

END SUBROUTINE OUTPUT_Surfaces


SUBROUTINE OUTPUT_Surface
    USE Global
    IMPLICIT NONE
    include "../tecio.F90"
    
    INTEGER:: I,J,K,L,M,N
    INTEGER:: Iav,Jav,Kav,count
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: NI2, NJ2, NK2


    INTEGER:: Ibgn,Iend,Jbgn,Jend,Kbgn,Kend
    INTEGER:: Itemp,Jtemp,Ktemp
    INTEGER:: surfibp,surfjbp,surfkbp
    INTEGER:: surfiep,surfjep,surfkep
    INTEGER:: surfip,surfjp,surfkp
    INTEGER:: Ind1,Ind2
    INTEGER:: IJK,MinorMax
    INTEGER:: NI_surf, NJ_surf, NK_surf
    INTEGER:: NI1_surf, NJ1_surf, NK1_surf
    
    INTEGER:: IJK_sur
    
    REAL:: P2,Cp,Cf,TT,RmiuW,RmiuW0,RmiuW1,Stan
    REAL:: Den,UU,VV,WW,VV0,VV1
    REAL:: Dn, TW,qw,ustar,uplus,yplus

    REAL:: Cp_av,Cf_av,T_av,Stan_av,PP_av
    REAL:: Den_av,UU_av,VV_av,WW_av
    REAL:: yplus_av
    
    INTEGER:: JJ,KK
    REAL:: TTve,KL,YY1,YY2,YY3,YY4,YY5,YY6
    REAL:: cve1,cve2,cve3,cve4,cve5,cve6
    REAL:: h1,h2,h3,h4,h5,h6          !各组分焓
    REAL:: Ds1,Ds2,Ds3,Ds4,Ds5,Ds6    !各组分质量扩散系数
    REAL:: cp1,cp2,cp3,cp4,cp5,cp6    !各组分定压比热容
    REAL:: KveL


    character(LEN=1):: NULLCHR

	character(LEN=100):: FL_surface
    character(LEN=100):: surfname
    CHARACTER(LEN=6):: cc
     
     Integer*4::   Debug,II,NIJK,VIsDouble
     INTEGER*4,ALLOCATABLE:: varlocation(:)
     INTEGER*4::   surf_VAR_num
      
      POINTER (NullPtr,Nulp)
      INTEGER*4:: Nulp(*)
      
      REAL*4,ALLOCATABLE:: VARs1(:,:,:)
      REAL*4,ALLOCATABLE:: VARs2(:,:,:)
      REAL*4,ALLOCATABLE:: VARs3(:,:,:)
      REAL*4,ALLOCATABLE:: VARs4(:,:,:)
      REAL*4,ALLOCATABLE:: VARs5(:,:,:)
      REAL*4,ALLOCATABLE:: VAR_value(:)

     NullPtr=0


    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    NI2=ThisBlock%NI2
    NJ2=ThisBlock%NJ2
    NK2=ThisBlock%NK2

!!!!!!
!	  output the flow infromation on the surface of PA
!!!!!!
!
	  do N=1,ThisBlock%Num_wall
	     Ibgn=ThisBlock%IJK_BE(N,1)
	     Iend=ThisBlock%IJK_BE(N,2)
	     Jbgn=ThisBlock%IJK_BE(N,3)
	     Jend=ThisBlock%IJK_BE(N,4)
	     Kbgn=ThisBlock%IJK_BE(N,5)
	     Kend=ThisBlock%IJK_BE(N,6)
!
         IJK_sur=ThisBlock%IJK_sur(N)
         FL_surface = ThisBlock%FL_sur (N)
         
         
         Ktemp=Kend-Kbgn+1
	     Jtemp=Jend-Jbgn+1
	     Itemp=Iend-Ibgn+1

         surfibp=0
         surfjbp=0
         surfkbp=0

         surfiep=1
         surfjep=1
         surfkep=1
         
         surfip =1
         surfjp =1
         surfkp =1	     
        
         if(IJK_sur == 1 ) then 
         surfiep=0
         surfip =0
         endif
       	 
         if(IJK_sur == 2 ) then 
         surfibp=1
         surfip =0
         endif
	     
         if(IJK_sur == 3 ) then
         surfjep=0
	     surfjp =0
         endif
         
         if(IJK_sur == 4 ) then
         surfjbp=1
         surfjp =0
         endif

	     if(IJK_sur == 5 ) then
         surfkep=0
	     surfkp =0
         endif
         
         if(IJK_sur == 6 ) then
         surfkbp=1
         surfkp =0
         endif
	


         NI_surf=Itemp+(surfiep-surfibp) 
         NJ_surf=Jtemp+(surfjep-surfjbp)
         NK_surf=Ktemp+(surfkep-surfkbp)

         NI1_surf=Itemp 
         NJ1_surf=Jtemp
         NK1_surf=Ktemp
!

      !WRITE(cc, '(I6.6)') KI

      
         WRITE(FL_surface, '(A,A)'), trim(FilePathPrefix), FL_surface
	         IF(KI >= KI_c)THEN
		        !DO Is=0,999
		          !IF( KI == INT(KI_c+Is*Ks) )THEN
		                !WRITE(cc, '(I3.3)') Is
		                WRITE(cc, '(I6.6)') (KI-KI_c)/Ksss

			            !WRITE(ThisBlock%FLPLT, '(A,A,A,A)') trim(FilePathPrefix), 'resu/RASST_', cc,".PLT"
                        WRITE(FL_surface, '(A,A,A,A)')trim(FL_surface),'_C', trim(cc), ".PLT"
		          !ENDIF
		        !ENDDO
	         ENDIF
       
        
         WRITE(*,*) 'output surf:', FL_surface
      
      
 	  !open(UNIT=17,FILE=trim(FilePathPrefix)//FL_surface//cc//'.PLT', MODE='WRITE',STATUS='UNKNOWN') !,SHARE='DENYWR')
         NULLCHR   = CHAR(0)
		 Debug     = 0						 ! 输出debug信息？0－不输出，1－输出少量，2－输出更多
		 VIsDouble = 0						 !是否采用双精度
!!!!!!
!	     以point方式输出
!!!!!!
!    >'X Y Z U V W  P T cp cf
!!!!!!
		surf_VAR_num =12
         ALLOCATE(varlocation(surf_VAR_num))
        
         varlocation(1:3)=1
		 varlocation(4:surf_VAR_num)=0	
         if(kind_field ==1)  varlocation(4:surf_VAR_num)=1		 



         I = TecIni100('SIMPLE DATASET'//NULLCHR, &
     &'X Y Z U V W P T Cp Cf Stan yplus'&
     &                  //NULLCHR, &
     &             trim(FL_surface)//NULLCHR, &
     &             '.'  //NULLCHR, &
     &             Debug, &
     &             VIsDouble)
		 I = TecZne100('Simple Zone'//NULLCHR,&
     &				0,&
     &              NI_surf,&
     &              NJ_surf,&
     &              NK_surf,&
     &				0,&
     &				0,&
     &				0,&
     &				1,&
     &				0,&
     &				0,&
     &				varlocation,&	
     &              nulp,&
     &                0)
   
    NIJK=NI_surf*NJ_surf*NK_surf
    
    ALLOCATE(VARs1(Iend+1,Jend+1,Kend+1))
    ALLOCATE(VARs2(Iend+1,Jend+1,Kend+1))
    ALLOCATE(VARs3(Iend+1,Jend+1,Kend+1))
    ALLOCATE(VARs4(Iend+1,Jend+1,Kend+1))
    ALLOCATE(VARs5(Iend+1,Jend+1,Kend+1))
    ALLOCATE(VAR_value(NIJK))

!!!!!!!!!!!	网格点，按照角点输出

     
    
    
      do K=Kbgn+surfkbp, Kend+surfkep
	  do J=Jbgn+surfjbp, Jend+surfjep
	  do I=Ibgn+surfibp, Iend+surfiep
          
          VARs1(I,J,K) = XX(I,J,K)
          VARs2(I,J,K) = YY(I,J,K)
          VARs3(I,J,K) = ZZ(I,J,K)
	  enddo
	  enddo
      enddo
!1
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn+surfibp:Iend+surfiep,Jbgn+surfjbp:Jend+surfjep,Kbgn+surfkbp:Kend+surfkep),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!2 
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn+surfibp:Iend+surfiep,Jbgn+surfjbp:Jend+surfjep,Kbgn+surfkbp:Kend+surfkep),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!3
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn+surfibp:Iend+surfiep,Jbgn+surfjbp:Jend+surfjep,Kbgn+surfkbp:Kend+surfkep),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!!!!!!!!!!!	流场变量，按照格心输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if(kind_field ==0) then
      NIJK=NI1_surf*NJ1_surf*NK1_surf

      do K=Kbgn,Kend
	  do J=Jbgn,Jend
	  do I=Ibgn,Iend
      	  Den    = V(1,I,J,K)
	      UU     = V(2,I,J,K)/Den
	      VV     = V(3,I,J,K)/Den
	      WW     = V(4,I,J,K)/Den
          VARs1(I,J,K) =UU
          VARs2(I,J,K) =VV
          VARs3(I,J,K) =WW
	  enddo
	  enddo
      enddo
!U
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!V
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!W
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)

     
      do K=Kbgn,Kend
	  do J=Jbgn,Jend
	  do I=Ibgn,Iend   
     
          VARs1(I,J,K) =PP(I,J,K)
          VARs2(I,J,K) =T(I,J,K)

      enddo
	  enddo
      enddo
!P
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!T
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)



      do K=Kbgn,Kend
	  do J=Jbgn,Jend
	  do I=Ibgn,Iend
      	  
          	P2 =PP(I,J,K)
            Cp =2.0*(P2-PPF)					!! Cp


	        Den = V(1,I,J,K)
            TT  =   T(I,J,K) 
            RmiuW =(1.+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
!            RmiuW =(1.+Csthlnd)/(TT+Csthlnd)*SQRT(TT*TT*TT)

	        UU=V(2,I,J,k)/V(1,I,J,k)
	        VV=V(3,I,J,k)/V(1,I,J,k)
	        WW=V(4,I,J,k)/V(1,I,J,k)
!
	   	    Dn = Dst(I,J,K)
            VV0=SQRT(V(2,I,J,k)*V(2,I,J,k)+   &
     &                 V(3,I,J,k)*V(3,I,J,k)+   &
     &                 V(4,I,J,k)*V(4,I,J,k))/V(1,I,J,K)
   
            Tw=RmiuW*VV0/DN/Ref
	        Cf=2.0_8*Tw
          
!            Cf= RmiuW 
           IF(V(2,I,J,K)<0.)Cf=-Cf

           ustar = sqrt(Tw/Den)
		   uplus = VV0/ustar
	       yplus = ustar*Dn/RmiuW*Den*Ref
          
           
 ! by luoj
  continue
  if (If_EquT) then
      
      TT  =   T(I,J,K)
  !    IF(IF_TwoTemp)then
  !    TT  =   T(I,J,K)
         
  !    ELSE
  !        qw=-(gam(I,J,K)/(gam(I,J,K)-1))*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
  !    ENDIF
  !    qw=-(gam(I,J,K)/(gam(I,J,K)-1))*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
      qw=-(1.4/(1.4-1))*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn     !by ydd
  else
     qw = 0.
  endif
     
     stan =qw/(1.4/0.4*PPF*(1.+0.2*XM*XM-Wall_Temp/Tinf))/Ref   !by ydd
  !   stan =qw/(gam(I,J,K)/(gam(I,J,K)-1)*PPF*(1.+0.5*(gam(I,J,K)-1)*XM*XM-Wall_Temp/Tinf))/Ref
    
    
          VARs1(I,J,K) =Cp
          VARs2(I,J,K) =Cf
          VARs3(I,J,K) =stan
          VARs4(I,J,K) =yplus

	  enddo
	  enddo
      enddo
!CP
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!Cf
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!St
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!yplus           
            var_value(1:NIJK)=RESHAPE(VARs4(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!格点输出!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else
      

      do K=Kbgn, Kend + surfkp
      do J=Jbgn, Jend + surfjp
	  do I=Ibgn, Iend + surfip
      	  
         UU_av=0
         VV_av=0
         WW_av=0
         count=0
         
         do kav = k-surfkp,k
	     do jav = j-surfjp,j
	     do iav = i-surfip,i
          
             Den    = V(1,Iav,Jav,Kav)
              if (Den > 0) then
	            UU     = V(2,Iav,Jav,Kav)/Den
	            VV     = V(3,Iav,Jav,Kav)/Den
	            WW     = V(4,Iav,Jav,Kav)/Den
          
                UU_av  = UU_av +UU
                VV_av  = VV_av +VV
                WW_av  = WW_av +WW
                count =  count +1
               endif
	  enddo
	  enddo
      enddo
          VARs1(I,J,K) =UU_av/count
          VARs2(I,J,K) =VV_av/count
          VARs3(I,J,K) =WW_av/count
      enddo
	  enddo
      enddo
!U
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!V
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!W
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)

     
      do K=Kbgn, Kend + surfkp
      do J=Jbgn, Jend + surfjp
	  do I=Ibgn, Iend + surfip  
     
         
         PP_av =0
         T_av  =0
         count =0
         
         do kav = k-surfkp,k
	     do jav = j-surfjp,j
	     do iav = i-surfip,i
          
             Den    = V(1,Iav,Jav,Kav)
              if (Den > 0) then
              PP_av =PP_av + PP(Iav,Jav,Kav)
              T_av  =T_av  + T (Iav,Jav,Kav)
              count =  count +1
          endif
          
          enddo
	      enddo
          enddo
              VARs1(I,J,K) =PP_av/count
              VARs2(I,J,K) =T_av /count
              
      
      enddo
	  enddo
      enddo
!P
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!T
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)



      do K=Kbgn, Kend + surfkp
      do J=Jbgn, Jend + surfjp
	  do I=Ibgn, Iend + surfip  
      	  
         Cp_av   =0
         Cf_av   =0
         stan_av =0
         yplus_av=0
         count   =0
         
         do kav = k-surfkp,k
	     do jav = j-surfjp,j
	     do iav = i-surfip,i
          
             Den    = V(1,Iav,Jav,Kav)
              if (Den > 0) then
            
            P2 =PP(Iav,Jav,Kav)
            Cp =2.0*(P2-PPF)					!! Cp

            TT  =   T(Iav,Jav,Kav)
            RmiuW =(1.+Csthlnd)/(TT*Csth1+Csthlnd)*(TT*Csth1)**1.5
!            RmiuW =(1.+Csthlnd)/(TT+Csthlnd)*SQRT(TT*TT*TT)

	        UU=V(2,Iav,Jav,Kav)/Den
	        VV=V(3,Iav,Jav,Kav)/Den
	        WW=V(4,Iav,Jav,Kav)/Den
!
	   	    Dn = Dst(Iav,Jav,Kav)
            VV0=SQRT(UU*UU+VV*VV+WW*WW)
            
            Tw=RmiuW*VV0/DN/Ref
	        Cf=2.0_8*Tw
            IF(V(2,I,J,K)<0.)Cf=-Cf

           ustar = sqrt(Tw/Den)
		   uplus = VV0/ustar
	       yplus = ustar*Dn/RmiuW*Den*Ref

  if (If_EquT) then
      
!      IF(IF_Chemical)then
      TT  =   T(Iav,Jav,Kav)
!      ELSE
!          qw=-(gam(Iav,Jav,Kav)/(gam(Iav,Jav,Kav)-1))*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
!      ENDIF
      qw=-(1.4_8/0.4_8)*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn     !by ydd
  
	   	        
     !qw=-(1.4_8/0.4_8)*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
!     qw=-(gam(Iav,Jav,Kav)/(gam(Iav,Jav,Kav)-1))*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
  else
     qw = 0.
  endif
   
        stan =qw/(1.4/0.4*PPF*(1.+0.2*XM*XM-Wall_Temp/Tinf))/Ref    !by ydd
!     stan =qw/(gam(Iav,Jav,Kav)/(gam(Iav,Jav,Kav)-1)*PPF*(1.+0.5*(gam(Iav,Jav,Kav)-1)*XM*XM-Wall_Temp/Tinf))/Ref
    
    
      Cp_av   =Cp_av    +Cp
      Cf_av   =Cf_av    +Cf
      stan_av =stan_av  +Stan
      yplus_av=yplus_av +yplus
      
      count = count+1
      endif
      
      enddo
	  enddo
      enddo
          
          VARs1(I,J,K) =Cp_av/count
          VARs2(I,J,K) =Cf_av/count
          VARs3(I,J,K) =stan_av/count
          VARs4(I,J,K) =yplus_av/count

	  enddo
	  enddo
      enddo
!CP
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!Cf
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!uplus
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!yplus           
            var_value(1:NIJK)=RESHAPE(VARs4(Ibgn:Iend+surfip,Jbgn:Jend+surfjp,Kbgn:Kend+surfkp),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)    
      
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif      
      
      
      
      
      I = TecEnd100()

      DEALLOCATE(VARs1)
      DEALLOCATE(VARs2)
      DEALLOCATE(VARs3)
      DEALLOCATE(VARs4)
      DEALLOCATE(VARs5)
      DEALLOCATE(var_value)
      DEALLOCATE(Varlocation)
!
	  enddo		!!DO for surfaces
!!!!!!


END SUBROUTINE OUTPUT_Surface
