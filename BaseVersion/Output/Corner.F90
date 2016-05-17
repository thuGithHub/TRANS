subroutine CornerPoints
use global
implicit none
integer::iblock
integer::I,J,K,NI1,NJ1,NK1,NI,NJ,NK,L

      do iblock=1,Max_Block
        call GLOBAL_SetPointers(iBlock)
        NI1=Thisblock%NI1
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1
        NI=NI1+1
        NJ=NJ1+1
        NK=NK1+1
       if(ThisBlock%ID_Present_blk==6)then
        open(unit=300,file="006.dat")
            do k=NK-2,NK+2
            do j=1,NJ-1
            do i=1,2
                write(300,*)KI,ThisBlock%ID_Present_blk,i,j,k,V(1,i,j,k),PP(i,j,k)
            enddo
            enddo
            enddo
        close(300)
        endif
       if(ThisBlock%ID_Present_blk==7)then
        open(unit=400,file="007.dat")
            do k=-2,3
            do j=1,NJ-1
            do i=1,2
                write(400,*)KI,ThisBlock%ID_Present_blk,i,j,k,V(1,i,j,k),PP(i,j,k)
            enddo
            enddo
            enddo
        close(400)
        endif

        do L=1,7
          do j=1,NJ1
                V(L,0,j,0)=1.5*V(L,0,j,1)-V(L,0,j,2)*0.5
                V(L,NI,j,0)=1.5*V(L,NI,j,1)-V(L,NI,j,2)*0.5
                V(L,0,j,NK)=1.5*V(L,0,j,NK1)-0.5*V(L,0,j,NK1-1)
                V(L,NI,j,NK)=1.5*V(L,NI,j,NK1)-0.5*V(L,NI,j,NK1-1)
          enddo
        enddo
        do L=1,9
          do j=1,NJ1
                dQdxyz(L,0,j,0)=1.5*dQdxyz(L,0,j,1)-dQdxyz(L,0,j,2)*0.5
                dQdxyz(L,NI,j,0)=1.5*dQdxyz(L,NI,j,1)-dQdxyz(L,NI,j,2)*0.5
                dQdxyz(L,0,j,NK)=1.5*dQdxyz(L,0,j,NK1)-0.5*dQdxyz(L,0,j,NK1-1)
                dQdxyz(L,NI,j,NK)=1.5*dQdxyz(L,NI,j,NK1)-0.5*dQdxyz(L,NI,j,NK1-1)
          enddo
        enddo

        do j=1,NJ1
                PP(0,j,0)=1.5*PP(0,j,1)-0.5*PP(0,j,2)
                PP(NI,j,0)=1.5*PP(NI,j,1)-0.5*PP(NI,j,2)
                PP(0,j,NK)=1.5*PP(0,j,NK1)-0.5*PP(0,j,NK1-1)
                PP(NI,j,NK)=1.5*PP(NI,j,NK1)-0.5*PP(NI,j,NK1-1)            
                T(0,j,0)=1.5*T(0,j,1)-0.5*T(0,j,2)
                T(NI,j,0)=1.5*T(NI,j,1)-0.5*T(NI,j,2)
                T(0,j,NK)=1.5*T(0,j,NK1)-0.5*T(0,j,NK1-1)
                T(NI,j,NK)=1.5*T(NI,j,NK1)-0.5*T(NI,j,NK1-1)            
                Rmiu(0,j,0)=1.5*Rmiu(0,j,1)-0.5*Rmiu(0,j,2)
                Rmiu(NI,j,0)=1.5*Rmiu(NI,j,1)-0.5*Rmiu(NI,j,2)
                Rmiu(0,j,NK)=1.5*Rmiu(0,j,NK1)-0.5*Rmiu(0,j,NK1-1)
                Rmiu(NI,j,NK)=1.5*Rmiu(NI,j,NK1)-0.5*Rmiu(NI,j,NK1-1)            
        enddo

        do L=1,7
          do k=0,NK
                V(L,0,0,k)=1.5*V(L,0,1,k)-0.5*V(L,0,2,k)
                V(L,0,NJ,k)=1.5*V(L,0,NJ1,k)-0.5*V(L,0,NJ1-1,k)
                V(L,NI,0,k)=1.5*V(L,NI,1,k)-0.5*V(L,NI,2,k)
                V(L,NI,NJ,k)=1.5*V(L,NI,NJ1,k)-0.5*V(L,NI,NJ1-1,k)
          enddo
        enddo
        do L=1,9
          do k=0,NK
                dQdxyz(L,0,0,k)=1.5*dQdxyz(L,0,1,k)-0.5*dQdxyz(L,0,2,k)
                dQdxyz(L,0,NJ,k)=1.5*dQdxyz(L,0,NJ1,k)-0.5*dQdxyz(L,0,NJ1-1,k)
                dQdxyz(L,NI,0,k)=1.5*dQdxyz(L,NI,1,k)-0.5*dQdxyz(L,NI,2,k)
                dQdxyz(L,NI,NJ,k)=1.5*dQdxyz(L,NI,NJ1,k)-0.5*dQdxyz(L,NI,NJ1-1,k)
          enddo
        enddo
        do k=0,NK
                PP(0,0,k)=1.5*PP(0,1,k)-0.5*PP(0,2,k)
                PP(0,NJ,k)=1.5*PP(0,NJ1,k)-0.5*PP(0,NJ1-1,k)
                PP(NI,0,k)=1.5*PP(NI,1,k)-0.5*PP(NI,2,k)
                PP(NI,NJ,k)=1.5*PP(NI,NJ1,k)-0.5*PP(NI,NJ1-1,k)        
                T(0,0,k)=1.5*T(0,1,k)-0.5*T(0,2,k)
                T(0,NJ,k)=1.5*T(0,NJ1,k)-0.5*T(0,NJ1-1,k)
                T(NI,0,k)=1.5*T(NI,1,k)-0.5*T(NI,2,k)
                T(NI,NJ,k)=1.5*T(NI,NJ1,k)-0.5*T(NI,NJ1-1,k)        
        
                Rmiu(0,0,k)=1.5*Rmiu(0,1,k)-0.5*Rmiu(0,2,k)
                Rmiu(0,NJ,k)=1.5*Rmiu(0,NJ1,k)-0.5*Rmiu(0,NJ1-1,k)
                Rmiu(NI,0,k)=1.5*Rmiu(NI,1,k)-0.5*Rmiu(NI,2,k)
                Rmiu(NI,NJ,k)=1.5*Rmiu(NI,NJ1,k)-0.5*Rmiu(NI,NJ1-1,k)        
        enddo
        
        do L=1,7
          do i=0,NI
            V(L,i,0,0)=1.5*V(L,i,1,0)-0.5*V(L,i,2,0)
            V(L,i,0,NK)=1.5*V(L,i,1,NK)-0.5*V(L,i,2,NK)
            V(L,i,NJ,0)=1.5*V(L,i,NJ1,0)-0.5*V(L,i,NJ1-1,0)
            V(L,i,NJ,NK)=1.5*V(L,i,NJ1,NK)-0.5*V(L,i,NJ1-1,NK)            
          enddo
        enddo
        do L=1,9
          do i=0,NI
            dQdxyz(L,i,0,0)=1.5*dQdxyz(L,i,1,0)-0.5*dQdxyz(L,i,2,0)
            dQdxyz(L,i,0,NK)=1.5*dQdxyz(L,i,1,NK)-0.5*dQdxyz(L,i,2,NK)
            dQdxyz(L,i,NJ,0)=1.5*dQdxyz(L,i,NJ1,0)-0.5*dQdxyz(L,i,NJ1-1,0)
            dQdxyz(L,i,NJ,NK)=1.5*dQdxyz(L,i,NJ1,NK)-0.5*dQdxyz(L,i,NJ1-1,NK)            
          enddo
        enddo

        do i=0,NI
            PP(i,0,0)=1.5*PP(i,1,0)-0.5*PP(i,2,0)
            PP(i,0,NK)=1.5*PP(i,1,NK)-0.5*PP(i,2,NK)
            PP(i,NJ,0)=1.5*PP(i,NJ1,0)-0.5*PP(i,NJ1-1,0)
            PP(i,NJ,NK)=1.5*PP(i,NJ1,NK)-0.5*PP(i,NJ1-1,NK)
            T(i,0,0)=1.5*T(i,1,0)-0.5*T(i,2,0)
            T(i,0,NK)=1.5*T(i,1,NK)-0.5*T(i,2,NK)
            T(i,NJ,0)=1.5*T(i,NJ1,0)-0.5*T(i,NJ1-1,0)
            T(i,NJ,NK)=1.5*T(i,NJ1,NK)-0.5*T(i,NJ1-1,NK)

            Rmiu(i,0,0)=1.5*Rmiu(i,1,0)-0.5*Rmiu(i,2,0)
            Rmiu(i,0,NK)=1.5*Rmiu(i,1,NK)-0.5*Rmiu(i,2,NK)
            Rmiu(i,NJ,0)=1.5*Rmiu(i,NJ1,0)-0.5*Rmiu(i,NJ1-1,0)
            Rmiu(i,NJ,NK)=1.5*Rmiu(i,NJ1,NK)-0.5*Rmiu(i,NJ1-1,NK)
        enddo
    enddo


end subroutine
        

