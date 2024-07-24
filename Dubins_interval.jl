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


function min_RS(orig_coord,orig_head,dest_coord,dest_head,r)    #   find the minimum path cost with max 2 segment path with RS structure
   
    # ds- only straight line segment
        
        theta1 = atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])     # heading angle of starting point

        theta1=zeroto2pi(theta1)    # converting the angle to zero to 2*pi

        theta2 = theta1     # for straight line segement start and ending heading is same
        
        if check_interval(theta1,theta2,orig_head,dest_head) # check if the heading angles lies withing the desired interval
        
            ds=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)    # cost of straight line segement
    
        else
            ds=1/0  # straight line segment doesn't exists within in given intervals the cost is set to infinity
        end
    
    ds=[ds,theta1,theta2]    # vector with cost and heading angles.
   
    
    
    #dr only curved segment

        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2) # distance between two targets

        dr=[1/0,0.0,0.0]  # cost is set to infinity

        if x_bar<=2*r   # curve segment can only exists for the case when distance between two targets is less than 2*r
        
            # now we would like to rotate the coordinates such that starting target is at (0,0) and ending target is at (x_bar,0)

            psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1]) # psi is the angle of rotation
            phi=asin(x_bar/(2*r)) # heading in rotated frame

            theta1=psi+phi  # rotating back to global frame
            theta2=psi-phi  # rotating back to global frame
            
            theta1=zeroto2pi(theta1)    # converting to 0 to 2*pi
            theta2=zeroto2pi(theta2)    # converting to 0 to 2*pi
            
            if check_interval(theta1,theta2,orig_head,dest_head) && dr[1]>=r*(zeroto2pi(theta1-theta2)) # checking the headings
            
                dr=r*(zeroto2pi(theta1-theta2)) # cost of path
                
                dr=[dr,theta1,theta2]

            end
            

        end
    
 
        # second solution for single circular segment
        if x_bar<=2*r
            
            phi=pi-asin(x_bar/(2*r))    # another possible solution for the circular arc
            psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1]) # psi is the angle of rotation
            
            theta1=psi+phi   # rotating back to global frame
            theta2=psi-phi   # rotating back to global frame
            
            theta1=zeroto2pi(theta1)    # converting to 0 to 2*pi
            theta2=zeroto2pi(theta2)    # converting to 0 to 2*pi

            if check_interval(theta1,theta2,orig_head,dest_head) && dr[1]>=r*(zeroto2pi(theta1-theta2))
            
                dr=r*(zeroto2pi(theta1-theta2)) # cost of path
                dr=[dr,theta1,theta2]    

            
            end

        end
        


    #RS1(theta1_min)  - 2 segment path starting with right hand arc from target 1 and then straight line segment to target 2 with heading of target 1 assumed to be theta1_min
        
    
        # now we would like to rotate the coordinates such that starting target is at (0,0) and ending target is at (x_bar,0)

        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)    # distance between two targets
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])      # find the angle of rotation to rotate the frame 
        
        
        theta1=orig_head[1] # heading of first target
        
        theta1_hat=theta1-psi   # heading of first target in rotated frame

        theta1_hat=zeroto2pi(theta1_hat)    # converted to 0 to 2pi

        solution=equation6(theta1_hat,x_bar,r)  # find the possible solutions based on the equation from the paper
     

        RS1_theta1min_value=1/0 # default cost is set to infinity
        RS1_theta1min=[1/0,0,0] # default cost is set to infinity and heading are set to 0

        for i in solution

     
            phi=i[1]    
            L=i[2]  # lenght of the straight line segment

            if L>=0 # to make sure straight line segment is positive

                theta2_hat=2*pi-phi # heading of ending target in rotated frame
        
                theta2_hat=zeroto2pi(theta2_hat)    # heading converted to 0 to 2pi

                theta2=theta2_hat+psi   # theta2 in rotated frame
        
                theta2=zeroto2pi(theta2)    # converting to 0 to 2pi

                 # Checking headings and make sure to store the minimum cost path among all possible solutions
                if check_interval(theta1,theta2,orig_head,dest_head) && RS1_theta1min_value>=r*zeroto2pi(theta1_hat+phi)+L 
                    RS1_theta1min_value=r*zeroto2pi(theta1_hat+phi)+L       # cost of the path
                    RS1_theta1min=[RS1_theta1min_value,theta1,theta2]       # cost and heading angles
            
                end

            end
  
        end
    

    # RS2(theta2_min) - 2 segment path starting with right hand arc from target 1 and then straight line segment to target 2 with heading of target 2 assumed to be theta2_min
        

        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)    # distance between two targets
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])      # find the angle of rotation to rotate the frame 
      
        theta2=dest_head[1] # heading of second target in global frame
        
        theta2_hat=theta2-psi   # heading of second target in rotated frame
        theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi

        phi=2*pi-theta2_hat # finding phi based on  heading of 2nd target in rotated frame
    
        phi=zeroto2pi(phi)  # converting 0 to 2pi
       
        RS1_theta2min_value=1/0 # setting default to  cost  to infinity
        RS1_theta2min=[1/0,0,0] # setting default to  cost to infinity and heading as zeroes

        if phi<=pi # phi should be less than equal to pi for solution to exists
        
            solution=equation7(phi,x_bar,r)  # find the possible solutions based on the equation from the paper

    
            for i in solution


                theta1_hat=i[1] # possible heading of target 1 in rotated frame
                L=i[2]      # length of straight line segment
                
                if L>=0 # length of straight line segment shoulde be greater than equal to zero
        
                    theta1=theta1_hat+psi       # heading of target 1 in global frame
        
                    theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
                    theta1=zeroto2pi(theta1)    # converting 0 to 2pi
        
                    if check_interval(theta1,theta2,orig_head,dest_head) && RS1_theta2min_value>=r*zeroto2pi(theta1_hat+phi)+L
        
                        RS1_theta2min_value=r*zeroto2pi(theta1_hat+phi)+L
                        RS1_theta2min=[RS1_theta2min_value,theta1,theta2]
            
                    end
 
                end
          
            end
        end

        
        # RS2(theta2_max)
            
        x_bar=sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)
        psi=atan(dest_coord[2] - orig_coord[2], dest_coord[1] - orig_coord[1])

        theta2=dest_head[2]

        theta2_hat=theta2-psi
        theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi

        phi=2*pi-theta2_hat

        phi=zeroto2pi(phi)    # converting 0 to 2pi

        RS1_theta2max_value=1/0
        RS1_theta2max=[1/0,0,0]

        if phi<=pi 

            solution=equation7(phi,x_bar,r)
    


            for i in solution
                theta1_hat=i[1]
                L=i[2]
            
                if L>=0 

                    theta1=theta1_hat+psi

                    theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
                    theta1=zeroto2pi(theta1)    # converting 0 to 2pi

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
    
        theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
    
        solution = equation1(theta1_hat,x_bar,r)
    
        RL1_theta1_min_value=1/0
        RL1_theta1_min=[1/0,0,0]
    
    
        for i in solution

            check=1
            
            if check==1
            theta2_hat = i[2]
            phi = i[1]
            
         
            theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
            phi=zeroto2pi(phi)    # converting 0 to 2pi
            theta2=theta2_hat+psi
    
            theta2=zeroto2pi(theta2)    # converting 0 to 2pi
    
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
    
        theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
    
        solution = equation1(theta1_hat,x_bar,r)

    
        RL1_theta1_max_value=1/0
        RL1_theta1_max=[1/0,0,0]
    
    
        for i in solution

            check=1
            
            if check==1
    
            theta2_hat = i[2]
            phi = i[1]
        
            theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
    
            phi=zeroto2pi(phi)    # converting 0 to 2pi
            theta2=theta2_hat+psi
    
            theta2=zeroto2pi(theta2)    # converting 0 to 2pi
    
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
     theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
    
     solution = equation2(theta2_hat,x_bar,r)

    
     RL2_theta2_min_value=1/0
     RL2_theta2_min=[1/0,0,0]
    
    
     for i in solution

        check=1
            
        if check==1
    
        theta1_hat = i[2]
        phi = i[1]
     
        theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
        phi=zeroto2pi(phi)    # converting 0 to 2pi
        theta1=theta1_hat+psi
    
        theta1=zeroto2pi(theta1)    # converting 0 to 2pi
    
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
     
     theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
    
     solution = equation2(theta2_hat,x_bar,r)
    
     RL2_theta2_max_value=1/0
     RL2_theta2_max=[1/0,0,0]
    
    
     for i in solution

        check=1
            
        if check==1
    
        theta1_hat = i[2]
        phi = i[1]
     
        theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
        phi=zeroto2pi(phi)    # converting 0 to 2pi
        theta1=theta1_hat+psi
    
        theta1=zeroto2pi(theta1)    # converting 0 to 2pi
    
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
    
                theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
                phi=zeroto2pi(phi)    # converting 0 to 2pi
                theta1_hat=theta2_hat
    
                theta1=theta1_hat+psi
    
                theta1=zeroto2pi(theta1)    # converting 0 to 2pi
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
    
            theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
            phi=zeroto2pi(phi)    # converting 0 to 2pi
            theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
        
            theta1=theta1_hat+psi
            theta1=zeroto2pi(theta1)    # converting 0 to 2pi
        
            theta2=theta2_hat+psi
            theta2=zeroto2pi(theta2)    # converting 0 to 2pi
    
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
    
            theta2_hat=zeroto2pi(theta2_hat)    # converting 0 to 2pi
            phi=zeroto2pi(phi)    # converting 0 to 2pi
            theta1_hat=zeroto2pi(theta1_hat)    # converting 0 to 2pi
    
            theta1=theta1_hat+psi
            theta1=zeroto2pi(theta1)    # converting 0 to 2pi
    
            theta2=theta2_hat+psi
            theta2=zeroto2pi(theta2)    # converting 0 to 2pi
    
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




function Dubins_interval(orig_coord,orig_head,dest_coord,dest_head,r)
 

    # 3 segment path 
    #Case1- When path is of the form RSR,RSL,LSR, LSL, RLR and LRL. 
    min_val = 1/0.0
    path_type = []

    for i in orig_head # iterate over extreme heading of first target

        for j in dest_head # iterate over extreme heading of second target
    
            temp1=copy(orig_coord)  # make a copy first target
            push!(temp1,i) # push the heading of first target
            temp2=copy(dest_coord)  # make a copy second target
            push!(temp2,j)  # push the heading of second target
            
            errcode, path = dubins_shortest_path(temp1,temp2, r)    # find the shortest path at the extreme headings
            
            #println(dubins_path_type(path))
            val1 = dubins_segment_length(path, 1)   # cost of first segment
            val2 = dubins_segment_length(path, 2)   # cost of second segment
            val3 = dubins_segment_length(path, 3)   # cost of third segment
            
                    
            if val1!=0.0 && val2!=0.0 && val3!=0 #  make sure to store 3 segment path only

                val=dubins_path_length(path)    

            else
                val=1/0.0
            end
            
            if val<min_val # we will store the minimum cost three segment path among all possible extreme headings

                global path1 = [val,i,j]    # storing cost and headings

            end
            
        end
    end
    
    
    
    
    if   sqrt((dest_coord[1]-orig_coord[1])^2+(dest_coord[2]-orig_coord[2])^2)<=4*r # 2 segment CC path is possible only when distance between two points is less than equalt to 4r
    
    
        #RL
        
            x1=min_RL_theta(orig_coord,orig_head,dest_coord,dest_head,r) # shortest path among, theta1=theta2, theta1+phi=pi or theta2+phi=pi
            x2=min_RL(orig_coord,orig_head,dest_coord,dest_head,r)      # shortest path either theta1 or theta2 is at exterme value
            

        #LR LR path has two ways to approach, either via reflection about x-axis or y-axis

        #LR1       reflection about x-axis
        
            temp1=copy(orig_coord)      #copy of starting point
            temp2=copy(dest_coord)      #copy of end point
            
            temp1[2]=-temp1[2]          #reflection about x-axis
            temp2[2]=-temp2[2]          #reflection about x-axis
            
            temp1_head=copy(orig_head)  # copy the heading interval
            temp2_head=copy(dest_head)  # copy the heading interval
            
            temp1_head[1]=zeroto2pi(-orig_head[2])  # new heading interval for reflected points
            temp1_head[2]=zeroto2pi(-orig_head[1])  # new heading interval for reflected points
            
            
            temp2_head[1]=zeroto2pi(-dest_head[2]) # new heading interval for reflected points
            temp2_head[2]=zeroto2pi(-dest_head[1])# new heading interval for reflected points
            
            x1_lr1=min_RL_theta(temp1,temp1_head,temp2,temp2_head,r) # shortest path among, theta1=theta2, theta1+phi=pi or theta2+phi=pi
            x2_lr1=min_RL(temp1,temp1_head,temp2,temp2_head,r)   # shortest path either theta1 or theta2 is at exterme value
            
            # now we will reflect the heading again to get back in global frame
            temp=x1_lr1[2]
            x1_lr1[2]=zeroto2pi(-x1_lr1[2]) 
            x1_lr1[3]=zeroto2pi(-x1_lr1[3])
            
            temp=x2_lr1[2]
            x2_lr1[2]=zeroto2pi(-x2_lr1[2])
            x2_lr1[3]=zeroto2pi(-x2_lr1[3])
            

        #LR2    reflection about y-axis  - reflection about x-axis or y-axis have same solution.
        
            temp1=copy(orig_coord) #copy of starting point
            temp2=copy(dest_coord) #copy of ending point
            
            temp1[1]=-temp1[1]   #reflection about y-axis
            temp2[1]=-temp2[1]   #reflection about y-axis 
            
            temp1_head=copy(orig_head) # copy thi heading
            temp2_head=copy(dest_head) # copy thi heading
            
            temp1_head[1]=zeroto2pi(-orig_head[2]) # reflection about y-axis
            temp1_head[2]=zeroto2pi(-orig_head[1]) # reflection about y-axis
            
            
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

    
    return x
    
end


    