clear;
clc;
close all

%%Declare Known Values
Mach_exit = 2.4;
P1= 1;
T1 = 298.15;
gamma = 1.4;

%Find nu at the Exit Mach number and the corresponding maximum theta value
nu_M = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(Mach_exit^2-1)))-atand(sqrt(Mach_exit^2-1));
theta_w_max = nu_M/2;

%Declare amount of characteristic lines and corresponding nodes
lines = 7;
nodes = .5*lines*(3+lines);

%Declare initial theta and delta_theta
theta_0 = .375;
delta_theta = 3;

%Find the values for each lines' first collision
for i = 1:lines
    theta(i) = theta_0 + (i-1)*delta_theta;
    nu(i) = theta(i);
    K_minus(i) = theta(i) + nu(i);
    K_plus(i) = theta(i) - nu(i);
end

%Assign crossover of values
next = lines+1;
theta(next) = theta(lines);
nu(next) = nu(lines);
K_minus(next) = K_minus(lines);
K_plus(next) = K_plus(lines);


%%For loop to find the values for the rest of the nodes

%Offset for matrix indices
start = 2;
ending = lines + 2;

for k=1:lines-1
    j=start;
    h=ending;
    
    %Assign values at the wall
    theta(h)=0;
    K_minus(h)=K_minus(j);
    nu(h)=K_minus(j)-theta(h);
    K_plus(h)=theta(h)-nu(h);

    %Calculate necessary values at each node where characteristic line
    %crosses other characteristic lines
    j=j+1;
    for i=h+1:lines-start+ending
        K_minus(i)=K_minus(j);
        K_plus(i)=K_plus(i-1);
        theta(i)=0.5*(K_plus(i)+K_minus(i));
        nu(i)=0.5*(K_minus(i)-K_plus(i));
        j=j+1;
    end
    if i==lines-start+ending
        h=i+1;
    else
        h=h+1;
    end
    
    %Assign values for characteristic line at top wall
    theta(h)=theta(h-1);
    nu(h)=nu(h-1);
    K_minus(h)=K_minus(h-1);
    K_plus(h)=K_plus(h-1);

    %Assign offset for next characteristic line
    start=start+1;
    ending=h+1;
end

%Determine Mach Number and mu value at each node based on reverse
%engineering Prandtl-Meyer equation
for i = 1:nodes
testM = 0;
testnu = -999;
while (abs(testnu-nu(i))) > 0.01
testnu= sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(testM^2-1)))-atand(sqrt(testM^2-1));
testM = testM + .0001;
end
M(i) = testM;
mu(i) = asind(1/M(i));
end

%Now that all values at each node have been found, we can plot the figure

%Find x and y values for node 1
x(1) = -1/tand(theta(1)-mu(1));
y(1) = 0;

%Find x and y-values for first characteristic line
for i = 2:lines
    x(i)=(1-y(i-1)+x(i-1)*tand(0.5*(mu(i-1)+theta(i-1)+mu(i)+theta(i))))/(tand(0.5*(mu(i-1)+theta(i-1)+mu(i)+theta(i)))-tand(theta(i)-mu(i)));
    y(i)=tand(theta(i)-mu(i))*x(i)+1;
end

%Find the first wall point
x(lines+1)=(y(lines)-1-x(lines)*tand((theta(lines)+theta(lines+1)+mu(lines)+mu(lines+1))*0.5))/(tand(0.5*(theta_w_max+theta(lines+1)))-tand((theta(lines)...
    +theta(lines+1)+mu(lines)+mu(lines+1))*0.5));
y(lines+1)=1+x(lines+1)*tand(0.5*(theta_w_max+theta(lines+1)));

i = 8;
k = lines;
%find x and y values for the rest of the nodes
for f = 2:lines
    for j = 1:k
        offset = j + i;
        %Centerline point
        if j == 1
            x(offset)=x(offset-k)-y(offset-k)/(tand(0.5*(theta(offset-k)+theta(offset)-mu(offset-k)-mu(offset))));
            y(offset)=0;
        %Wall Point
        else if j == k
        x(offset)=(x(offset-k)*tand(0.5*(theta(offset-k)+theta(offset)))-y(offset-k)+y(offset-1)-x(offset-1)*tand((theta(offset-1)+theta(offset)+mu(offset-1)+mu(offset))*0.5))/(tand(0.5*(theta(offset-k)...
            +theta(offset)))-tand((theta(offset-1)+theta(offset)+mu(offset-1)+mu(offset))*0.5));
        y(offset)=y(offset-k)+(x(offset)-x(offset-k))*tand(0.5*(theta(offset-k)+theta(offset))); 
        %Points in between
        else
            slope_last = tand(0.5*(theta(offset)+theta(offset-1)+mu(offset)+mu(offset-1)));
            slope_across = tand(0.5*(theta(offset)+theta(offset-k)-mu(offset)-mu(offset-k)));
            x(offset)=(y(offset-k)-y(offset-1)+slope_last*x(offset-1)-slope_across*x(offset-k))/(slope_last-slope_across);
            y(offset)=y(offset-1)+(x(offset)-x(offset-1))*slope_last;
        end
        end
    end
    i = offset;
    k = k-1;     
end

m = lines + 1;
offset = m-1;
%Plot the figure, beginning with the wall
for i = 1:lines
    if i == 1
        plot([0 x(m)], [1, y(m)], 'color', 'black')
        hold on
    else
        plot([x(m-offset-1),x(m)], [y(m-offset-1),y(m)], 'color', 'black')
        hold on
    end
    m = m + offset;
    offset = offset - 1;
end

%Plot the initial lines from the inlet edge to the centerline for each
%characteristic line
for i = 1:lines
    value = i;
    subtraction = 0;
    curr = value;
    while value > 0
        if value == i
        plot([0,x(value)], [1,y(value)], 'color', 'black')
        hold on    
        else
        plot([x(last), x(curr)], [y(last), y(curr)], 'color', 'black')    
        end
        last = curr;
        curr = last + (lines - subtraction);
        value = value-1;
        subtraction = subtraction + 1;
    end
end
first = 1;
%Plot the rebounded lines throught each node
for i = 1:lines
    value = lines - i;
    for j = first:(first + value)
        plot([x(j), x(j+1)], [y(j), y(j+1)], 'color', 'black')
    end
    first = first + value + 2;
end


%Plot characteristics
xlim([0, x(nodes) + 1])
ylim([0, y(nodes) + 1])
xlabel("Dimensionless x value")
ylabel("Dimensionless y value")
title("Minimum Length Nozzle for Exit Mach Number of 2.4")

%Label each node
for i = 1:nodes
    txt = [num2str(i)];
    text(x(i)-0.03, y(i)+ 0.075, txt);
end

%Create node matrix for table labels
for i = 1:nodes
node(i) = i;
end

%Transpose all matrices to be read in table
node = transpose(node);
M = transpose(M);
K_minus = transpose(K_minus);
K_plus = transpose(K_plus);
theta = transpose(theta);
nu = transpose(nu);
mu = transpose(mu);

%Create and display table
table = table(node, K_minus, K_plus, theta, nu, M, mu);
display(table)

