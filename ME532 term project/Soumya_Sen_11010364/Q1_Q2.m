clc; clear all;

%-------------------------------------------------------------------
%FEM solution
%-------------------------------------------------------------------


%----------------------------------------------
%Inputs
%----------------------------------------------
syms x;
syms xi;
%Loading type - Distributed force f(X)
prompt = 'Input coefficients of f(x) = a0 + a1x + a2x^2 + a3x^3: '
prompt = 'a0 ';
a0 = input(prompt);
prompt = 'a1 ';
a1 = input(prompt);
prompt = 'a2 ';
a2 = input(prompt);
prompt = 'a3 ';
a3 = input(prompt);

f = a0 + a1*x + a2*x^2 + a3*x^3 

%Material property - EA(x)
prompt = 'Input coefficients of EA(x) = b0 + b1x + b2x^2: '
prompt = 'b0 ';
b0 = input(prompt);
prompt = 'b1 ';
b1 = input(prompt);
prompt = 'b2 ';
b2 = input(prompt);

EA = b0 + b1*x + b2*x^2

%Number of elements - n_e
prompt = 'Input the number of elements: ';
n_e = input(prompt);

%Length of the rod - l
prompt = 'Input the length of the rod: ';
l = input(prompt);

%Boundary Conditions - bc

check1 = 0; check2 = 0; check3 = 0; check4 = 0; %Initializing parameters to check inputs given at boundary

%At x = 0
prompt = 'At x = 0, Enter 0 for Dirichlet, 1 for Neuman: ';
bc1 = input(prompt);
if bc1 == 0
	prompt = 'Enter displacement at x = 0: ';
	u_1 = input(prompt);
    check1 = 1;
elseif bc1 == 1
	prompt = 'Enter force at x = 0: ';
	P_1 = input(prompt);
    check2 = 1;
end

%At x = L
prompt = 'At x = L, Enter 0 for Dirichlet, 1 for Neuman: ';
bcL = input(prompt);
if bcL == 0
	prompt = 'Enter displacement at x = L: ';
	u_L = input(prompt);
    check3 = 1;
elseif bcL == 1
	prompt = 'Enter force at x = L: ';
	P_L = input(prompt);
    check4 = 1;
end


%Lagrangian polynomial - p
prompt = 'Order of Lagrangian polynomial: ';
p = input(prompt);

count1=0;
K_g(n_e*p+1,n_e*p+1) = 0;
F_g(n_e*p+1,1) = 0;


%Meshing
l_e = l/n_e; %Length of each element

%NINT
%NINT = ceil((p+1)/2);
NINT = 5;
if NINT == 1
    fid = fopen('gauss_points1.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 2
    fid = fopen('gauss_points2.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 3
    fid = fopen('gauss_points3.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 4
    fid = fopen('gauss_points4.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 5
    fid = fopen('gauss_points5.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 6
    fid = fopen('gauss_points6.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 7
    fid = fopen('gauss_points7.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 8
    fid = fopen('gauss_points8.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 9
    fid = fopen('gauss_points9.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
elseif NINT == 10
    fid = fopen('gauss_points10.txt', 'r');
    gauss_points = fscanf(fid, '%f %f\n', [2 NINT])';
    fclose(fid);
else
    error('Change p');
end


cord = zeros(n_e,2);
for index1 = 1:n_e
    
    product = 1;
    
	cord(index1, 1) = count1;    %Storing the cordinates of the elements
	cord(index1, 2) = count1+l_e;
	count1 = cord(index1, 2);


	%--------------------------------------------
	%Processor
	%--------------------------------------------

	xi_j = zeros(p+1,1);
	for index2 = 1:p+1
		%Locating the points for pth order polynomial - xi_j
		xi_j(index2,1) = -1+((2/p)*(index2-1));
	end


	%x -> xi
	x = ((1-xi)/2)*cord(index1,1) + ((1+xi)/2)*cord(index1,2);
	EA = b0 + b1*x + b2*x^2;
    f = a0 + a1*x + a2*x^2 + a3*x^3;

	%N = zeros(p+1,1);

	%Shape functions and derivatives
	for index3 = 1:p+1
        product = 1;
		for index4 = 1:p+1

			if index3 ~= index4;
				N(index3,1) = product*((xi - xi_j(index4,1))/(xi_j(index3,1) - xi_j(index4,1)));
				product = N(index3,1);
            end
            
        end 
        dN(index3,1) = diff(N(index3,1),xi);
    end
    
    K_e = zeros(p+1, p+1); f_e = zeros(p+1, 1);
	%Element stiffness matrix
	for index5 = 1:p+1
		for index6 = 1:p+1
            
			for q = 1:NINT
                w_i = gauss_points(q, 1); zeta = gauss_points(q, 2);
                Ni_1 = subs( dN(index5,1), xi , zeta);
                Ni_2 = subs( dN(index6,1), xi , zeta);
                EA_1 = subs( EA, xi , zeta);
                K_e(index5,index6) = K_e(index5,index6) + ((2/l_e)* EA_1* (Ni_1) *(Ni_2)*w_i);
			end
		end
		for q = 1:NINT
            w_i = gauss_points(q, 1); zeta = gauss_points(q, 2);
            f_1 = subs( f, xi , zeta);
			f_e(index5,1) = f_e(index5,1) + ((l_e/2)*f_1*(subs( N(index5,1), xi , zeta)*w_i));
		end
    end
    
    %Assembly
    K_g((index1-1)*p + 1:index1*p + 1,(index1-1)*p + 1:index1*p + 1) = K_g((index1-1)*p + 1:index1*p + 1,(index1-1)*p + 1:index1*p + 1) + K_e(1:p+1,1:p+1);
    F_g((index1-1)*p + 1:index1*p + 1,1) = F_g((index1-1)*p + 1:index1*p + 1,1) + f_e(1:p+1,1);

end


%Applying Boundary conditions
if bc1 == 0
	%Dirichlet condition
	K_g(1,:) = 0;
	K_g(:,1) = 0;
	K_g(1,1) = 1;
	F_g(1) = u_1;
elseif bc1 == 1
	%Neuman condition
	F_g(1) = F_g(1) + P_1;
end

if bcL == 0
	%Dirichlet condition
	K_g(n_e*p+1,:) = 0;
	K_g(:,n_e*p+1) = 0;
	K_g(n_e*p+1,n_e*p+1) = 1;
	F_g(n_e*p+1,1) = u_L;
elseif bcL == 1
	%Neuman condition
	F_g(n_e*p+1,1) = F_g(n_e*p+1,1) + P_L;
end

%Finding the FEM solution
u_g = (K_g)\F_g

%The strain energy
s_e = 0.5*u_g'*K_g*u_g


%-------------------------------------------------------------------
%Actual solution
%-------------------------------------------------------------------

syms u(x);
EA = b0 + b1*x + b2*x^2;
f = a0 + a1*x + a2*x^2 + a3*x^3;

if (check1 == 1) && (check4 == 1)
u_act = dsolve(diff(EA*diff(u)) + f == 0, u(0) == u_1, subs( EA*diff(u), x , l) == P_L)
elseif (check1 == 1) && (check3 == 1)
u_act = dsolve(diff(EA*diff(u)) + f == 0, u(0) == u_1, u(l) == u_L)
elseif (check2 == 1) && (check3 == 1)
u_act = dsolve(diff(EA*diff(u)) + f == 0, subs( EA*diff(u), x , 0) == P_1, u(l) == u_L)
elseif (check2 == 1) && (check4 == 1)
u_act = dsolve(diff(EA*diff(u)) + f == 0, subs( EA*diff(u), x , 0) == P_1, subs( EA*diff(u), x , l) == P_L)
end


%---------------------------------------------------------------------
%Plotting the values
%---------------------------------------------------------------------

figure
ezplot(u_act,[0,l/(n_e*p),l])
hold on;
 
x = 0:l/(n_e*p):l;
plot(x, u_g, 'ro')

xlabel('x (nodes)');      
ylabel('u (displacement)');
title('Plot of u-fem and u-actual');
legend('u-act','u-fem');

hold off