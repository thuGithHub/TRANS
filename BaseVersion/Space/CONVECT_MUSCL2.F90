SUBROUTINE CONVECT_MUSCL2

    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    
    REAL:: UU,VV,WW,PP0,TT
    REAL:: DR,DL
    real::rrr

    INTERFACE
        REAL FUNCTION MUSCL2(var1,var2)  !muscl2
                    REAL:: var1,var2
        END FUNCTION
    END INTERFACE
    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1



!        VL(3,ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm         
       DO i=0,NI+1
       DO j=0,NJ+1
	   DO k=0,NK+1
	   DO L=1,5
		  VR(1,L,i,j,k)=VL(1,L,i,j,k)
		  VR(2,L,i,j,k)=VL(2,L,i,j,k)
		  VR(3,L,i,j,k)=VL(3,L,i,j,k)
	   end do
	   end do
	   end do
       end do


!        i direction
         Do k=1,NK1
         Do j=1,NJ1
         Do i=0,NI

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(1,1,i+1,j,k)-UU                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=UU-VR(1,1,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(1,1,i+1,j,k)=UU+MUSCL2(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+MUSCL2(dr,dl)

	   VR(1,1,i,j,k)=UU-MUSCL2(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(1,2,i+1,j,k)-VV                 ! dr=vil(i+1,j,k)-  vy(i,j,k)
         DL=VV-VR(1,2,i,j,k)                   ! dl=vy(i,j,k)    -vir(i,j,k) 

         VL(1,2,i+1,j,k)=VV+MUSCL2(DR,DL)    ! vil(i+1,j,k)=vy(i,j,k)+MUSCL2(dr,dl)

	   VR(1,2,i,j,k)=VV-MUSCL2(DL,DR)      ! vir(i,j,k)  =vy(i,j,k)-MUSCL2(dl,dr)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(1,3,i+1,j,k)-WW                 ! dr=(wil(i+1,j,k)-  vz(i,j,k))
         DL=WW-VR(1,3,i,j,k)                   ! dl=(vz(i,j,k)    -wir(i,j,k))

         VL(1,3,i+1,j,k)=WW+MUSCL2(DR,DL)    ! wil(i+1,j,k)=vz(i,j,k)+MUSCL2(dr,dl)

	   VR(1,3,i,j,k)=WW-MUSCL2(DL,DR)      ! wir(i,j,k)  =vz(i,j,k)-MUSCL2(dl,dr)



!         limite P


         PP0=PP(i,j,k)

         DR=(VL(1,4,i+1,j,k)-PP0) !*2.*PP(i,j,k)/(VL(1,4,i+1,j,k)+PP0)                 ! dr=(pil(i+1,j,k)-  p(i,j,k))
         DL=(PP0-VR(1,4,i,j,k))   !*2.*PP(i,j,k)/(PP0+VR(1,4,i,j,k))                   ! dl=(p(i,j,k)    -pir(i,j,k))

         VL(1,4,i+1,j,k)=PP0+MUSCL2(DR,DL)    ! pil(i+1,j,k)=p(i,j,k)+MUSCL2(dr,dl)

	   VR(1,4,i,j,k)=PP0-MUSCL2(DL,DR)      ! pir(i,j,k)  =p(i,j,k)-MUSCL2(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=(VL(1,5,i+1,j,k)-TT)!*2.*TT/(VL(1,5,i+1,j,k)+TT)                 ! dr=(til(i+1,j,k)-  tp(i,j,k))
         DL=(TT-VR(1,5,i,j,k))  !*2.*TT/(TT+VR(1,5,i,j,k))                   ! dl=(tp(i,j,k)    -tir(i,j,k))

         VL(1,5,i+1,j,k)=TT+MUSCL2(DR,DL)    ! til(i+1,j,k)=tp(i,j,k)+MUSCL2(dr,dl)

	   VR(1,5,i,j,k)=TT-MUSCL2(DL,DR)      ! tir(i,j,k)  =tp(i,j,k)-MUSCL2(dl,dr)
        
        !limit rrc, by ydd
        rrr=rad(i,j,k)
        DR=rccl(1,i+1,j,k)-rrr
        DL=rrr-rccr(1,i,j,k)

        rccl(1,i+1,j,k)=rrr+MUSCL2(DR,DL)
        rccr(1,i,j,k)=rrr-MUSCL2(DL,DR)

	   end do
	   end do
	   end do



!        j direction

	   Do k=1,NK1
         Do j=0,NJ
         Do i=1,NI1

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(2,1,i,j+1,k)-UU                 !  dr=ujl(i,j+1,k)-  vx(i,j,k)
         DL=UU-VR(2,1,i,j,k)                   !  dl=vx(i,j,k)    -ujr(i,j,k))

         VL(2,1,i,j+1,k)=UU+MUSCL2(DR,DL)    ! ujl(i,j+1,k)=vx(i,j,k)+MUSCL2(dr,dl)

	   VR(2,1,i,j,k)=UU-MUSCL2(DL,DR)      ! ujr(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(2,2,i,j+1,k)-VV                 !  dr=(vjl(i,j+1,k)-  vy(i,j,k))
         DL=VV-VR(2,2,i,j,k)                   !  dl=(vy(i,j,k)    -vjr(i,j,k))

         VL(2,2,i,j+1,k)=VV+MUSCL2(DR,DL)    ! vjl(i,j+1,k)=vy(i,j,k)+MUSCL2(dr,dl)

	   VR(2,2,i,j,k)=VV-MUSCL2(DL,DR)      ! vjr(i,j,k)  =vy(i,j,k)-MUSCL2(dl,dr)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(2,3,i,j+1,k)-WW                 ! dr=(wjl(i,j+1,k)-  vz(i,j,k))
         DL=WW-VR(2,3,i,j,k)                   ! dl=(vz(i,j,k)    -wjr(i,j,k))

         VL(2,3,i,j+1,k)=WW+MUSCL2(DR,DL)    ! wjl(i,j+1,k)=vz(i,j,k)+MUSCL2(dr,dl)

	   VR(2,3,i,j,k)=WW-MUSCL2(DL,DR)      ! wjr(i,j,k)  =vz(i,j,k)-MUSCL2(dl,dr)

    
!         limite P


         PP0=PP(i,j,k)

         DR=(VL(2,4,i,j+1,k)-PP0)!*2.*PP0/(VL(2,4,i,j+1,k)+PP0)                 ! dr=(pjl(i,j+1,k)-  p(i,j,k))
         DL=(PP0-VR(2,4,i,j,k))  !*2.*PP0/(VR(2,4,i,j,k)+PP0)                   ! dl=(p(i,j,k)    -pjr(i,j,k))

         VL(2,4,i,j+1,k)=PP0+MUSCL2(DR,DL)    ! pjl(i,j+1,k)=p(i,j,k)+MUSCL2(dr,dl)

	   VR(2,4,i,j,k)=PP0-MUSCL2(DL,DR)      ! pjr(i,j,k)  =p(i,j,k)-MUSCL2(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=(VL(2,5,i,j+1,k)-TT) !*2.*TT/(VL(2,5,i,j+1,k)+TT)                 ! dr=(tjl(i,j+1,k)-  tp(i,j,k))
         DL=(TT-VR(2,5,i,j,k))   !*2.*TT/(VR(2,5,i,j,k)+TT)                   ! dl=(tp(i,j,k)    -tjr(i,j,k))

         VL(2,5,i,j+1,k)=TT+MUSCL2(DR,DL)    ! tjr(i,j,k)  =tp(i,j,k)-MUSCL2(dl,dr)

	   VR(2,5,i,j,k)=TT-MUSCL2(DL,DR)      ! tjl(i,j+1,k)=tp(i,j,k)+MUSCL2(dr,dl)

        !limit rrc, by ydd
        rrr=rad(i,j,k)
        DR=rccl(2,i,j+1,k)-rrr
        DL=rrr-rccr(2,i,j,k)

        rccl(2,i,j+1,k)=rrr+MUSCL2(DR,DL)
        rccr(2,i,j,k)=rrr-MUSCL2(DL,DR)

	   end do
	   end do
	   end do


!        k direction
         Do k=0,NK
         Do j=1,NJ1
         Do i=1,NI1

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(3,1,i,j,k+1)-UU                 ! dr=ukl(i,j,k+1)-  vx(i,j,k)
         DL=UU-VR(3,1,i,j,k)                   ! dl=vx(i,j,k)    -ukr(i,j,k)

         VL(3,1,i,j,k+1)=UU+MUSCL2(DR,DL)    ! ukr(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)

	   VR(3,1,i,j,k)=UU-MUSCL2(DL,DR)      ! ukl(i,j,k+1)=vx(i,j,k)+MUSCL2(dr,dl)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(3,2,i,j,k+1)-VV                 ! dr=(vkl(i,j,k+1)-  vy(i,j,k))
         DL=VV-VR(3,2,i,j,k)                   ! dl=(vy(i,j,k)    -vkr(i,j,k))

         VL(3,2,i,j,k+1)=VV+MUSCL2(DR,DL)    ! vkr(i,j,k)  =vy(i,j,k)-MUSCL2(dl,dr)

	   VR(3,2,i,j,k)=VV-MUSCL2(DL,DR)      ! vkl(i,j,k+1)=vy(i,j,k)+MUSCL2(dr,dl)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(3,3,i,j,k+1)-WW                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=WW-VR(3,3,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,3,i,j,k+1)=WW+MUSCL2(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+MUSCL2(dr,dl)

	   VR(3,3,i,j,k)=WW-MUSCL2(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)



!         limite P


         PP0=PP(i,j,k)

         DR=(VL(3,4,i,j,k+1)-PP0) !*2.*PP0/(VL(3,4,i,j,k+1)+PP0)                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=(PP0-VR(3,4,i,j,k))   !*2.*PP0/(VR(3,4,i,j,k)+PP0)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,4,i,j,k+1)=PP0+MUSCL2(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+MUSCL2(dr,dl)

	   VR(3,4,i,j,k)=PP0-MUSCL2(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=(VL(3,5,i,j,k+1)-TT)!*2.*TT/(VL(3,5,i,j,k+1)+TT)                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=(TT-VR(3,5,i,j,k))  !*2.*TT/(TT+VR(3,5,i,j,k))                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,5,i,j,k+1)=TT+MUSCL2(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+MUSCL2(dr,dl)

	     VR(3,5,i,j,k)=TT-MUSCL2(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-MUSCL2(dl,dr)

        !limit rrc, by ydd
        rrr=rad(i,j,k)
        DR=rccl(3,i,j,k+1)-rrr
        DL=rrr-rccr(3,i,j,k)

        rccl(3,i,j,k+1)=rrr+MUSCL2(DR,DL)
        rccr(3,i,j,k)=rrr-MUSCL2(DL,DR)

	   end do
	   end do
	   end do



END SUBROUTINE CONVECT_MUSCL2
