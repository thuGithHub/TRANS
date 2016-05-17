SUBROUTINE Average_Surfaces
    USE Global
    IMPLICIT NONE
    
    INTEGER:: iBlock
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        CALL Average_Surface
    ENDDO

END SUBROUTINE Average_Surfaces


SUBROUTINE Average_Surface
    USE Global
    IMPLICIT NONE
    include "../tecio.F90"
    
    INTEGER:: I,J,K,L,M,N
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: NI2, NJ2, NK2
    
    INTEGER:: NI_surf, NJ_surf, NK_surf
    INTEGER:: NI1_surf, NJ1_surf, NK1_surf

    INTEGER:: Ibgn,Iend,Jbgn,Jend,Kbgn,Kend
    INTEGER:: Itemp,Jtemp,Ktemp
    INTEGER:: surfibp,surfjbp,surfkbp
    INTEGER:: surfiep,surfjep,surfkep
    INTEGER:: IJK,MinorMax
    
    INTEGER:: IJK_sur
      
    INTEGER::KZ
    REAL::   TIME

    REAL:: P2,Cp,Cf,TT,RmiuW,RmiuW0,RmiuW1,ustar,uplusb,yplusb,Cfb,Stan
    REAL:: Den,UU,VV,WW,VV0,VV1
    REAL:: Dn, TW,qw
    REAL:: Vab_04,Vab_05,Vab_06,Vab_07,Vab_08,Vab_09,Prms

	character(LEN=1):: NULLCHR
    character(LEN=100):: FL_surface
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
!!
		 IJK_sur=ThisBlock%IJK_sur(N)
         FL_surface = ThisBlock%FL_suraver(N)
        

         
	     Ktemp=Kend-Kbgn+1
	     Jtemp=Jend-Jbgn+1
	     Itemp=Iend-Ibgn+1

         surfibp=0
         surfjbp=0
         surfkbp=0

         surfiep=1
         surfjep=1
         surfkep=1
	     
        
         if(IJK_sur == 1 ) surfiep=0
       	 if(IJK_sur == 2 ) surfibp=1
   
	     if(IJK_sur == 3 ) surfjep=0
	     if(IJK_sur == 4 ) surfjbp=1


	     if(IJK_sur == 5 ) surfkep=0
	     if(IJK_sur == 6 ) surfkbp=1
	


         NI_surf=Itemp+(surfiep-surfibp) 
         NJ_surf=Jtemp+(surfjep-surfjbp)
         NK_surf=Ktemp+(surfkep-surfkbp)

         NI1_surf=Itemp 
         NJ1_surf=Jtemp
         NK1_surf=Ktemp
!

 !
         
         KZ=MOD(KI,KQ)
	     IF(Kz == 0 )THEN
 	   
!
	  

!	     以point方式输出
!!!!!!
!    >'X Y Z U V W  P T cp cf Stan yplus prms
!!!!!!
		surf_VAR_num =13
         ALLOCATE(varlocation(surf_VAR_num))
        
         varlocation(1:3)=1
		 varlocation(4:surf_VAR_num)=0	
         
         NULLCHR   = CHAR(0)
		 Debug     = 0						 ! 输出debug信息？0－不输出，1－输出少量，2－输出更多
		 VIsDouble = 0						 !是否采用双精度
!!!!!!	 



       I = TecIni100('SIMPLE DATASET'//NULLCHR, &
     &'X Y Z U V W P T Cp Cf Stan yplus prms'&
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
!!!!!!!!!!!	流场变量，按照格心输出
	  
      NIJK=NI1_surf*NJ1_surf*NK1_surf

      do K=Kbgn,Kend
	  do J=Jbgn,Jend
	  do I=Ibgn,Iend

          VARs1(I,J,K) =Ub(I,J,K)
          VARs2(I,J,K) =Vb(I,J,K)
          VARs3(I,J,K) =Wb(I,J,K)
          VARs4(I,J,K) =pb(I,J,K)
          VARs5(I,J,K) =tb(I,J,K)
	  enddo
	  enddo
      enddo
!4
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!5
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!6
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)

!7
            var_value(1:NIJK)=RESHAPE(VARs4(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!8
            var_value(1:NIJK)=RESHAPE(VARs5(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)

      do K=Kbgn,Kend
	  do J=Jbgn,Jend
	  do I=Ibgn,Iend
      	  
      
          
         	TT  = Tb(I,J,K)
            RmiuW =(1.+Csthlnd)/(TT*Csth1+Csthlnd)*(TT*Csth1)**1.5
	        Den = Rb(I,J,K)	        
			
			UU=Ub(I,J,K)
	        VV=Vb(I,J,K)
	        WW=Wb(I,J,K)
            Dn = Dst(I,J,K)
            VV0=SQRT(UU*UU+VV*VV+WW*WW)
   
            Tw=RmiuW*VV0/DN/Ref
	        Cfb=2.0_8*Tw
              IF(UU<0.)Cfb=-Cfb

           ustar = sqrt(Tw/Den)
	       uplusb = VV0/ustar
	       yplusb = ustar*Dn/RmiuW*Den*Ref

  if (If_EquT) then
	   	        
     qw=-(1.4_8/0.4_8)*PPF*(RmiuW/0.72)*(TT-Wall_Temp/Tinf)/Dn
  else
     qw = 0.
  endif
   
     stan =qw/(1.4/0.4*PPF*(1.+0.2*XM*XM-Wall_Temp/Tinf))/Ref
 
         
          VARs1(I,J,K) = 2.0* (pb(I,J,K)-PPF)  !cp
          VARs2(I,J,K) = Cfb                   !cf 
          VARs3(I,J,K) = stan                 !stan
          VARs4(I,J,K) = yplusb               !yplus
          VARs5(I,J,K) = pprms(I,J,K)


	  enddo
	  enddo
      enddo
!9
			var_value(1:NIJK)=RESHAPE(VARs1(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!10
			var_value(1:NIJK)=RESHAPE(VARs2(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!11
			var_value(1:NIJK)=RESHAPE(VARs3(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!12           
            var_value(1:NIJK)=RESHAPE(VARs4(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)
!13
            var_value(1:NIJK)=RESHAPE(VARs5(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend),(/NIJK/))
			II      = TecDat100(NIJK,var_value,VIsDouble)


      I = TecEnd100()
         
      DEALLOCATE(VARs1)
      DEALLOCATE(VARs2)
      DEALLOCATE(VARs3)
      DEALLOCATE(VARs4)
      DEALLOCATE(VARs5)
      DEALLOCATE(var_value)
      DEALLOCATE(Varlocation)

      
      endif
      
      enddo

END SUBROUTINE Average_Surface




