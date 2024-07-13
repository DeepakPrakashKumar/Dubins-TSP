using Dubins


function check_interval(angle1,angle2,orig_head,dest_head) # to check if the heading at starting and destination lies within the allowed intervals

    
    if orig_head[2]>=orig_head[1]       

        if angle1>=orig_head[1] && orig_head[2]>=angle1 
            a1=1.0
        else
            a1=0.0
        end
    else
        
        if (angle1>=orig_head[1] && angle1<=2*pi) || (angle1>=0.0 && angle1<=orig_head[2])
            a1=1.0
        else
            a1=0.0
        end
    end

    if dest_head[2]>=dest_head[1]

        if angle2>=dest_head[1] && dest_head[2]>=angle2 
            a2=1.0
        else
            a2=0.0
        end
    else
        
        if (angle2>=dest_head[1] && angle2<=2*pi) || (angle2>=0.0 && angle2<=dest_head[2])
            a2=1.0
        else
            a2=0.0
        end
    end


    return a2==1.0 && a1==1

end



function zeroto2pi(angle)       # converts angle between 0 to 2 pi
    
    if angle >= 2*pi
        angle -= 2 * π
    end
    
    if angle < 0
        angle += 2 * π
    end

    return angle
end


function equation1(theta1,x_bar,r)      # equations to find theta2 and phi for RL at extreme values of theta1
        
    sol=[]

    a=4*r^2*cos(theta1)
    b=4*r^2*sin(theta1)-4*x_bar*r
    c=4*r^2+x_bar^2-2*x_bar*r*sin(theta1)

    beta=atan(b,a)
    if abs(c/sqrt(a^2+b^2))<=1
        phi=acos(c/sqrt(a^2+b^2))-beta
        theta2=atan(x_bar/r-2*sin(phi)-sin(theta1),2*cos(phi)-cos(theta1))
        push!(sol,[phi theta2])     # first possible solution

        phi=2*pi-acos(c/sqrt(a^2+b^2))-beta
        theta2=atan(x_bar/r-2*sin(phi)-sin(theta1),2*cos(phi)-cos(theta1))
        push!(sol,[phi theta2])     # second possible solution

    end
    return sol



end

function equation2(theta2,x_bar,r)      # equations to find theta1 and phi for RL at extreme values of theta2

    sol=[]

    a=4*r^2*cos(theta2)
    b=4*r^2*sin(theta2)-4*x_bar*r
    c=4*r^2+x_bar^2-2*x_bar*r*sin(theta2)

    beta=atan(b,a)
    if abs(c/sqrt(a^2+b^2))<=1
        phi=acos(c/sqrt(a^2+b^2))-beta
        theta1=atan(x_bar/r-2*sin(phi)-sin(theta2),2*cos(phi)-cos(theta2))
        push!(sol,[phi theta1])     # first possible solution

        phi=2*pi-acos(c/sqrt(a^2+b^2))-beta
        theta1=atan(x_bar/r-2*sin(phi)-sin(theta2),2*cos(phi)-cos(theta2))
        push!(sol,[phi theta1])     # second possible solution

    end
    return sol
end

function equation3(x_bar,r)  # equations to find theta1, theta2 and phi  for RL when theta1=theta2
        
    sol=[]



    if abs(x_bar/(4*r))<=1

            theta1=asin(x_bar/(4*r))
            phi=atan(x_bar/(2*r)-sin(theta1),cos(theta1))
            push!(sol,[phi, theta1, theta1])        # first possible solution

            theta1=pi-asin(x_bar/(4*r))
            phi=atan(x_bar/(2*r)-sin(theta1),cos(theta1))
            push!(sol,[phi, theta1, theta1])        # second possible solution

    end

    return sol

end


function equation4(x_bar,r) # equations to find theta1, theta2 and phi  for RL when theta1+phi=pi
    
    sol=[]

    if abs((x_bar^2+8*r^2)/(6*x_bar*r))<=1

        theta1=asin((x_bar^2+8*r^2)/(6*x_bar*r))
        phi=pi-theta1
        theta2=atan(x_bar/r-3*sin(theta1),-3*cos(theta1))
        push!(sol,[phi, theta1, theta2])        # first possible solution

        theta1=pi-asin((x_bar^2+8*r^2)/(6*x_bar*r))
        phi=pi-theta1
        theta2=atan(x_bar/r-3*sin(theta1),-3*cos(theta1))
        push!(sol,[phi, theta1, theta2])    # second possible solution

     end    


    return sol

end

function equation5(x_bar,r)  # equations to find theta1, theta2 and phi  for RL when theta2+phi=pi

    sol=[]

    if abs((x_bar^2+8*r^2)/(6*x_bar*r))<=1

        theta2=asin((x_bar^2+8*r^2)/(6*x_bar*r))
        phi=pi-theta2
        theta1=atan(x_bar/r-3*sin(theta2),-3*cos(theta2))
        push!(sol,[phi, theta1, theta2])        # first possible solution


        theta2=pi-asin((x_bar^2+8*r^2)/(6*x_bar*r))
        phi=pi-theta2
        theta1=atan(x_bar/r-3*sin(theta2),-3*cos(theta2))
        push!(sol,[phi, theta1, theta2])        # second possible solution


     end   

    return sol
end



function equation6(theta1,x_bar,r) # equations to find phi and L  for RS at theta1_min
 
    sol=[]

    if x_bar^2-2*x_bar*r*sin(theta1)>=0
        
        L=sqrt(x_bar^2-2*x_bar*r*sin(theta1))
        beta=acos(r/sqrt(r^2+L^2))

        phi=zeroto2pi(acos((r/sqrt(r^2+L^2))*cos(theta1))-beta)
        
        if phi<=pi
            push!(sol,[phi,L])
        end
       
        phi=zeroto2pi(2*pi-acos((r/sqrt(r^2+L^2))*cos(theta1))-beta)
       
        if phi<=pi
            push!(sol,[phi,L])      # first solution
        end

    end

    return sol

end

function equation7(phi,x_bar,r)     # equations to find theta1 and L  for RS at at extreme values of theta2
 
    sol=[]

    a=1
    b=2*r*sin(phi)*cos(phi)-2*x_bar*cos(phi)-2*r*sin(phi)*cos(phi)
    c=x_bar^2-2*x_bar*r*sin(phi)

    if b^2-4*a*c>=0
        L1=(-b-sqrt(b^2-4*a*c))/(2*a)
        L2=(-b+sqrt(b^2-4*a*c))/(2*a)

        if L1>=0
            theta1=atan(-sin(phi)-L1*cos(phi)/r+x_bar/r,cos(phi)-L1*sin(phi)/r)
            push!(sol,[theta1,L1])      #first solution
        end

        
        if L2>=0
            theta1=atan(-sin(phi)-L2*cos(phi)/r+x_bar/r,cos(phi)-L2*sin(phi)/r)
            push!(sol,[theta1,L2])      #second solution
        end

    end


    return sol

end


function min_RS(orig_coord,orig_head,dest_coord,dest_head,r)
   
    # ds
        
        theta1 = atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])

        theta1=zeroto2pi(theta1)

        theta2 = theta1
        
        if check_interval(theta1,theta2,orig_head,dest_head) 
        
        ds=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
    
        else
            ds=1/0
        end
    
    ds=[ds,theta1,theta2]    
   
    
    
    #dr
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        dr=[1/0,theta1,theta2]
        if x_bar<=2*r
        
        phi=asin(x_bar/(2*r))
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
        
        theta1=psi+phi
        theta2=psi-phi
        
        theta1=zeroto2pi(theta1)
        theta2=zeroto2pi(theta2)
        if check_interval(theta1,theta2,orig_head,dest_head) && dr[1]>=r*(zeroto2pi(theta1-theta2))
        
        dr=r*(zeroto2pi(theta1-theta2))
            
        dr=[dr,theta1,theta2]

        end
        

        end
    
 

    if x_bar<=2*r
        
        phi=pi-asin(x_bar/(2*r))
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
        
        theta1=psi+phi
        theta2=psi-phi
        
        theta1=zeroto2pi(theta1)
        theta2=zeroto2pi(theta2)

        if check_interval(theta1,theta2,orig_head,dest_head) && dr[1]>=r*(zeroto2pi(theta1-theta2))
        
            dr=r*(zeroto2pi(theta1-theta2))
            dr=[dr,theta1,theta2]    

        
        end

    end
    


    #RS1(theta1_min)
        
         x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
         psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
        
        
        theta1=orig_head[1]
        
        theta1_hat=theta1-psi

        theta1_hat=zeroto2pi(theta1_hat)

        solution=equation6(theta1_hat,x_bar,r)
     

        RS1_theta1min_value=1/0
        RS1_theta1min=[1/0,0,0]

        for i in solution

            check=1
            
            if check==1

            phi=i[1]
            L=i[2]
            if L>=0 

                theta2_hat=2*pi-phi
        
                theta2_hat=zeroto2pi(theta2_hat)

                theta2=theta2_hat+psi
        
                theta2=zeroto2pi(theta2)
        
                if check_interval(theta1,theta2,orig_head,dest_head) && RS1_theta1min_value>=r*zeroto2pi(theta1_hat+phi)+L
        
                    RS1_theta1min_value=r*zeroto2pi(theta1_hat+phi)+L
                    RS1_theta1min=[RS1_theta1min_value,theta1,theta2]
            
                end

            end
        end
        end
    

    # RS2(theta2_min)
        
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
        
        theta2=dest_head[1]
        
        theta2_hat=theta2-psi
        theta2_hat=zeroto2pi(theta2_hat)

        phi=2*pi-theta2_hat
    
        phi=zeroto2pi(phi)
       
        RS1_theta2min_value=1/0
        RS1_theta2min=[1/0,0,0]

        if phi<=pi 
        
            solution=equation7(phi,x_bar,r)


    
            for i in solution

                check=1
                if check==1

                theta1_hat=i[1]
                L=i[2]
                
                if L>=0 
        
                    theta1=theta1_hat+psi
        
                    theta1_hat=zeroto2pi(theta1_hat)
                    theta1=zeroto2pi(theta1)
        
                    if check_interval(theta1,theta2,orig_head,dest_head) && RS1_theta2min_value>=r*zeroto2pi(theta1_hat+phi)+L
        
                        RS1_theta2min_value=r*zeroto2pi(theta1_hat+phi)+L
                        RS1_theta2min=[RS1_theta2min_value,theta1,theta2]
            
                    end
 
                end
                end
            end
        end

        
        # RS2(theta2_max)
            
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])

        theta2=dest_head[2]

        theta2_hat=theta2-psi
        theta2_hat=zeroto2pi(theta2_hat)

        phi=2*pi-theta2_hat

        phi=zeroto2pi(phi)

        RS1_theta2max_value=1/0
        RS1_theta2max=[1/0,0,0]

        if phi<=pi 

            solution=equation7(phi,x_bar,r)
    


            for i in solution
                theta1_hat=i[1]
                L=i[2]
            
                if L>=0 

                    theta1=theta1_hat+psi

                    theta1_hat=zeroto2pi(theta1_hat)
                    theta1=zeroto2pi(theta1)

                    if check_interval(theta1,theta2,orig_head,dest_head) && RS1_theta2max_value>=r*zeroto2pi(theta1_hat+phi)+L

                        RS1_theta2max_value=r*zeroto2pi(theta1_hat+phi)+L
                        RS1_theta2max=[RS1_theta2max_value,theta1,theta2]
        
                    end

                end
            end

        end



            R=[ds,dr,RS1_theta1min,RS1_theta2min,RS1_theta2max]
        
            min_index = argmin(v -> v[1], R)
            
        return min_index
            
        
        
end


    
    
    
function min_RL(orig_coord,orig_head,dest_coord,dest_head,r)
        
        # RL1 theta1 min
    
     
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
       
        theta1=orig_head[1]
       
        theta1_hat=theta1-psi
    
        theta1_hat=zeroto2pi(theta1_hat)
    
        solution = equation1(theta1_hat,x_bar,r)
    
        RL1_theta1_min_value=1/0
        RL1_theta1_min=[1/0,0,0]
    
    
        for i in solution

            check=1
            
            if check==1
            theta2_hat = i[2]
            phi = i[1]
            
         
            theta2_hat=zeroto2pi(theta2_hat)
            phi=zeroto2pi(phi)
            theta2=theta2_hat+psi
    
            theta2=zeroto2pi(theta2)
    
            if check_interval(theta1,theta2,orig_head,dest_head) && RL1_theta1_min_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
        
                RL1_theta1_min_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
                RL1_theta1_min=[RL1_theta1_min_value,theta1,theta2]
                
            end
        
        end
        
        end
        
        # RL1 theta1 max
     
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
       
        theta1=orig_head[2]
       
        theta1_hat=theta1-psi
    
        theta1_hat=zeroto2pi(theta1_hat)
    
        solution = equation1(theta1_hat,x_bar,r)

    
        RL1_theta1_max_value=1/0
        RL1_theta1_max=[1/0,0,0]
    
    
        for i in solution

            check=1
            
            if check==1
    
            theta2_hat = i[2]
            phi = i[1]
        
            theta2_hat=zeroto2pi(theta2_hat)
    
            phi=zeroto2pi(phi)
            theta2=theta2_hat+psi
    
            theta2=zeroto2pi(theta2)
    
            if check_interval(theta1,theta2,orig_head,dest_head)  && RL1_theta1_max_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
        
               RL1_theta1_max_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
               RL1_theta1_max=[RL1_theta1_max_value,theta1,theta2]
        
    
            end
    
        end  
        end
    
     # RL2 theta1 min
    
     initial_guess = [0.0, 0.0]
    
     
     x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
     psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
    
     theta2=dest_head[1]
    
     theta2_hat=theta2-psi
     theta2_hat=zeroto2pi(theta2_hat)
    
     solution = equation2(theta2_hat,x_bar,r)

    
     RL2_theta2_min_value=1/0
     RL2_theta2_min=[1/0,0,0]
    
    
     for i in solution

        check=1
            
        if check==1
    
        theta1_hat = i[2]
        phi = i[1]
     
        theta1_hat=zeroto2pi(theta1_hat)
        phi=zeroto2pi(phi)
        theta1=theta1_hat+psi
    
        theta1=zeroto2pi(theta1)
    
        if check_interval(theta1,theta2,orig_head,dest_head)  && RL2_theta2_min_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
     
                RL2_theta2_min_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
                RL2_theta2_min=[RL2_theta2_min_value,theta1,theta2]
        end
    end
    end
    
     # RL2 theta2 max
    
     initial_guess = [0.0, 0.0]
    
    
     x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
     psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
    
     theta2=dest_head[2]
    
     theta2_hat=theta2-psi
     
     theta2_hat=zeroto2pi(theta2_hat)
    
     solution = equation2(theta2_hat,x_bar,r)
    
     RL2_theta2_max_value=1/0
     RL2_theta2_max=[1/0,0,0]
    
    
     for i in solution

        check=1
            
        if check==1
    
        theta1_hat = i[2]
        phi = i[1]
     
        theta1_hat=zeroto2pi(theta1_hat)
        phi=zeroto2pi(phi)
        theta1=theta1_hat+psi
    
        theta1=zeroto2pi(theta1)
    
        if check_interval(theta1,theta2,orig_head,dest_head)  && RL2_theta2_max_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
     
            RL2_theta2_max_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
            RL2_theta2_max=[RL2_theta2_max_value,theta1,theta2]
            
        end
    end
    end
    
     R=[RL1_theta1_min,RL1_theta1_max,RL2_theta2_min,RL2_theta2_max]
    
     min_index = argmin(v -> v[1], R)
     
         return min_index
        
end



    
function min_RL_theta(orig_coord,orig_head,dest_coord,dest_head,r)
    
        # theta1=theta2
        RL_theta1=[1/0,0,0]
        RL_theta1_phi=[1/0,0,0]
        RL_theta2_phi=[1/0,0,0]
    
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])
           
        solution=equation3(x_bar,r)

        RL_theta1_value=1/0
        RL_theta1=[1/0,0,0]
        
        for i in solution   
            check=1

            if check==1

                theta2_hat=i[3]
                phi=i[1]
    
                theta2_hat=zeroto2pi(theta2_hat)
                phi=zeroto2pi(phi)
                theta1_hat=theta2_hat
    
                theta1=theta1_hat+psi
    
                theta1=zeroto2pi(theta1)
                theta2=theta1
    
                RL_theta1_value=1/0
                RL_theta1=[1/0,0,0]
    
                if check_interval(theta1,theta2,orig_head,dest_head) && RL_theta1_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
     
                    RL_theta1_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
                    RL_theta1=[RL_theta1_value,theta1,theta2]
                
        
                end
            end 
        end
        # theta1+phi=pi 
        
        solution=equation4(x_bar,r)

        
        RL_theta1_phi_value=1/0
        RL_theta1_phi=[1/0,0,0]

        for i in solution

            check=1

            if check==1

        
            theta2_hat = i[3]
            phi = i[1]
            theta1_hat = i[2]
    
            theta2_hat=zeroto2pi(theta2_hat)
            phi=zeroto2pi(phi)
            theta1_hat=zeroto2pi(theta1_hat)
        
            theta1=theta1_hat+psi
            theta1=zeroto2pi(theta1)
        
            theta2=theta2_hat+psi
            theta2=zeroto2pi(theta2)
    
            RL_theta1_phi_value=1/0
            RL_theta1_phi=[1/0,0,0]
    
            if check_interval(theta1,theta2,orig_head,dest_head) && RL_theta1_phi_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
     
                RL_theta1_phi_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
                RL_theta1_phi=[RL_theta1_value,theta1,theta2]
            end
        end
        end
    
        # theta2+phi=pi 
    
        solution=equation5(x_bar,r)
        RL_theta2_phi_value=1/0
        RL_theta2_phi=[1/0,0,0]
            
        for i in solution
            check=1

            if check==1

            theta2_hat = i[3]
            phi = i[1]
            theta1_hat = i[2]
    
            theta2_hat=zeroto2pi(theta2_hat)
            phi=zeroto2pi(phi)
            theta1_hat=zeroto2pi(theta1_hat)
    
            theta1=theta1_hat+psi
            theta1=zeroto2pi(theta1)
    
            theta2=theta2_hat+psi
            theta2=zeroto2pi(theta2)
    
            RL_theta2_phi_value=1/0
            RL_theta2_phi=[1/0,0,0]
    
            if check_interval(theta1,theta2,orig_head,dest_head) && RL_theta2_phi_value>=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
     
                RL_theta2_phi_value=r*(zeroto2pi(theta1_hat+phi)+zeroto2pi(phi+theta2_hat))
                RL_theta2_phi=[RL_theta1_value,theta1,theta2]
    
            end
        end
        end 
    
    R=[RL_theta1,RL_theta1_phi,RL_theta2_phi]
   
    min_index = argmin(v -> v[1], R)
    
        return min_index
end


#Case1- When path is of the form RSR,RSL,LSR, LSL, RLR and LRL.

function Dubins_interval(orig_coord,orig_head,dest_coord,dest_head,r)
 

    # 3 segment path 
    min_val = 1/0.0
    path_type = []

    for i in orig_head

        for j in dest_head
    
            temp1=copy(orig_coord)
            push!(temp1,i)
            temp2=copy(dest_coord)
            push!(temp2,j)
            
            errcode, path = dubins_shortest_path(temp1,temp2, r)
            
            #println(dubins_path_type(path))
            val1 = dubins_segment_length(path, 1)
            val2 = dubins_segment_length(path, 2)
            val3 = dubins_segment_length(path, 3)
            
                    
            if val1!=0.0 && val2!=0.0 && val3!=0

                val=dubins_path_length(path)

            else
                val=1/0.0
            end
            
            if val<min_val

                min_val = val
                path_type = dubins_path_type(path)
                global path1 = [val,i,j]
            end
            
        end
    end
    
    
    
    
    if   sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)<=4*r # 2 segment CC path is possible only when distance between two points is less than equalt to 4r
    
    
        #RL
        
            x1=min_RL_theta(orig_coord,orig_head,dest_coord,dest_head,r) # shortest path among, theta1=theta2, theta1+phi=pi or theta2+phi=pi
            x2=min_RL(orig_coord,orig_head,dest_coord,dest_head,r)      # shortest path either theta1 or theta2 is at exterme value
            

        #LR

        #LR1       reflection about x-axis
        
            temp1=copy(orig_coord)      #copy of starting point
            temp2=copy(dest_coord)      #copy of end point
            
            temp1[2]=-temp1[2]          #reflection about x-axis
            temp2[2]=-temp2[2]          #reflection about x-axis
            
            temp1_head=copy(orig_head)
            temp2_head=copy(dest_head)
            
            temp1_head[1]=zeroto2pi(-orig_head[2])  # new heading interval for reflected points
            temp1_head[2]=zeroto2pi(-orig_head[1])
            
            
            temp2_head[1]=zeroto2pi(-dest_head[2])
            temp2_head[2]=zeroto2pi(-dest_head[1])
            
            x1_lr1=min_RL_theta(temp1,temp1_head,temp2,temp2_head,r)
            x2_lr1=min_RL(temp1,temp1_head,temp2,temp2_head,r)
            
            temp=x1_lr1[2]
            x1_lr1[2]=zeroto2pi(-x1_lr1[2])
            x1_lr1[3]=zeroto2pi(-x1_lr1[3])
            
            temp=x2_lr1[2]
            x2_lr1[2]=zeroto2pi(-x2_lr1[2])
            x2_lr1[3]=zeroto2pi(-x2_lr1[3])
            
        #LR2    reflection about y-axis  - reflection about x-axis or y-axis have same solution.
        
            temp1=copy(orig_coord)
            temp2=copy(dest_coord)
            
            temp1[1]=-temp1[1]
            temp2[1]=-temp2[1]
            
            temp1_head=copy(orig_head)
            temp2_head=copy(dest_head)
            
            temp1_head[1]=zeroto2pi(-orig_head[2])
            temp1_head[2]=zeroto2pi(-orig_head[1])
            
            
            temp2_head[1]=zeroto2pi(-dest_head[2])
            temp2_head[2]=zeroto2pi(-dest_head[1])
            
            x1_lr2=min_RL_theta(temp2,temp2_head,temp1,temp1_head,r)
            x2_lr2=min_RL(temp2,temp2_head,temp1,temp1_head,r)


            temp=x1_lr2[2]
            x1_lr2[2]=zeroto2pi(-x1_lr2[3])
            x1_lr2[3]=zeroto2pi(-temp)
            
            temp=x2_lr2[2]
            x2_lr2[2]=zeroto2pi(-x2_lr2[3])
            x2_lr2[3]=zeroto2pi(-temp)
            
            
        
    else

        x1=[1/0,0,0]
        x2=[1/0,0,0]
        x1_lr1=[1/0,0,0]
        x1_lr2=[1/0,0,0]
        x2_lr1=[1/0,0,0]
        x2_lr2=[1/0,0,0]
    
    end
    
    #RS    
    
         x3=min_RS(orig_coord,orig_head,dest_coord,dest_head,r) #shortest RS path 
    
    #SR     relection about y-axis
    
        temp1=copy(orig_coord)
        temp2=copy(dest_coord)
        
        temp1[1]=-temp1[1]
        temp2[1]=-temp2[1]
        
        temp1_head=copy(orig_head)
        temp2_head=copy(dest_head)
        
        temp1_head[1]=zeroto2pi(-orig_head[2])
        temp1_head[2]=zeroto2pi(-orig_head[1])
        
        
        temp2_head[1]=zeroto2pi(-dest_head[2])
        temp2_head[2]=zeroto2pi(-dest_head[1])
        
        x3_sr=min_RS(temp2,temp2_head,temp1,temp1_head,r)
        
        
        temp=x3_sr[2]
        x3_sr[2]=zeroto2pi(-x3_sr[3])
        x3_sr[3]=zeroto2pi(-temp)
        
        
    #LS     reflection about x-axis
    
        temp1=copy(orig_coord)
        temp2=copy(dest_coord)
        
        temp1[2]=-temp1[2]
        temp2[2]=-temp2[2]
        
        temp1_head=copy(orig_head)
        temp2_head=copy(dest_head)
        
        temp1_head[1]=zeroto2pi(-orig_head[2])
        temp1_head[2]=zeroto2pi(-orig_head[1])
        
        
        temp2_head[1]=zeroto2pi(-dest_head[2])
        temp2_head[2]=zeroto2pi(-dest_head[1])
        
        x3_ls=min_RS(temp1,temp1_head,temp2,temp2_head,r)
        
    

        x3_ls[2]=zeroto2pi(-x3_ls[2])
        x3_ls[3]=zeroto2pi(-x3_ls[3])
        

    
    #SL reflection about x-axis then y-axis

        temp1=copy(orig_coord)
        temp2=copy(dest_coord)
        
        temp1[2]=-temp1[2]
        temp2[2]=-temp2[2]
        
        temp1[1]=-temp1[1]
        temp2[1]=-temp2[1]
        
        temp1_head=copy(orig_head)
        temp2_head=copy(dest_head)
        
        
        x3_sl=min_RS(temp2,temp2_head,temp1,temp1_head,r)



        temp=x3_sl[2]
        x3_sl[2]=zeroto2pi(x3_sl[3])
        x3_sl[3]=zeroto2pi(temp)
        
        
    
    
    x=[x1,x2,x3,path1,x1_lr1,x1_lr2,x2_lr1,x2_lr2,x3_sr,x3_ls,x3_sl]
    
    minimum = argmin(v -> v[1], x)      # minimum among all possible paths

    
    return minimum
    
end


    