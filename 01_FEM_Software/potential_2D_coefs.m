% b: semi-chord
% a: elastic axis position
% c: hinge line

function amat = potential_2D_coefs(rho,U,b,a,c)


% matrices related to [h, alpha, beta]' and first and second derivatives
A1_nc(1,1) = -(pi*rho*b^2);
A1_nc(1,2) = -(-pi*rho*b^2*b*a);
A1_nc(2,1) = pi*rho*b^2*b*a;
A1_nc(2,2) = -pi*rho*b^2*b^2*(1/8+a^2);

A1_nc(1,3) = -(-rho*b^3*(c*acos(c)-1/3*(2+c^2)*sqrt(1-c^2)));
A1_nc(2,3) = -rho*b^4*((1/8+c^2)*acos(c) - 1/8*c*sqrt(1-c^2)*(7+2*c^2)+...
    (c-a)*(1/3*sqrt(1-c^2)*(2+c^2)-c*acos(c)));
A1_nc(3,1) = rho*b^3*(c*acos(c)-1/3*(2+c^2)*sqrt(1-c^2));
A1_nc(3,2) = -rho*b^4*((1/8+c^2)*acos(c)-1/8*c*sqrt(1-c^2)*(7+2*c^2)+...
    (c-a)*(1/3*(2+c^2)*sqrt(1-c^2)-c*acos(c)));
A1_nc(3,3) = rho*b^4/pi*(1/4*c*sqrt(1-c^2)*acos(c)*(7+2*c^2)-...
    (1/8+c^2)*(acos(c))^2 -1/8*(1-c^2)*(5*c^2+4));

A1_c(1,1) = 0;
A1_c(1,2) = 0;
A1_c(2,1) = 0;
A1_c(2,2) = 0;

A1_c(1,3) = 0;
A1_c(2,3) = 0;
A1_c(3,1) = 0;
A1_c(3,2) = 0;
A1_c(3,3) = 0;

A2_nc(1,1) = 0;    
A2_nc(1,2) = -(pi*rho*b^2*U);
A2_nc(2,1) = pi*rho*b^2*U;
A2_nc(2,2) = 0;

A2_nc(1,3) = -(-rho*b^2*U*(c*sqrt(1-c^2)-acos(c)));
A2_nc(2,3) = -rho*b^3*U*(1/3*sqrt(1-c^2)*(c^2-1)-(c-a)*(c*sqrt(1-c^2)-acos(c)));
A2_nc(3,1) = -rho*b^2*U*(c*sqrt(1-c^2)-acos(c));
A2_nc(3,2) = rho*b^3*U*(a*(c*sqrt(1-c^2)-acos(c))+1/3*(sqrt(1-c^2))^3-...
    1/3*(2+c^2)*sqrt(1-c^2)+c*acos(c));
     
A2_c(1,1) = -(2*pi*rho*U*b*1);
A2_c(1,2) = -(2*pi*rho*U*b*b*(1/2-a));
A2_c(2,1) = -2*pi*rho*U*b^2*1/2*1 + b*(a+1/2)*2*pi*rho*U*b*1;
A2_c(2,2) = -2*pi*rho*U*b^2*1/2*b*(1/2-a) + b*(a+1/2)*2*pi*rho*U*b*b*(1/2-a);

A2_c(1,3) = -(2*pi*rho*U*b*b/(2*pi)*((1-2*c)*acos(c)+(2-c)*sqrt(1-c^2)));
A2_c(2,3) = -2*pi*rho*U*b^2*(1/2-(a+1/2))*b/(2*pi)*((1-2*c)*acos(c)+(2-c)*sqrt(1-c^2));
fat1 = -2*rho*U*b^2*(1/2*(acos(c)-c*sqrt(1-c^2)) + ...
    ((1+c/2)*sqrt(1-c^2)-(c+1/2)*acos(c))); % common factor 1
A2_c(3,1) = fat1;
A2_c(3,2) = fat1*b*(1/2-a);
A2_c(3,3) = fat1*b/(2*pi)*((1-2*c)*acos(c)+(2-c)*sqrt(1-c^2));

A3_nc(1,1) = 0;
A3_nc(1,2) = 0;
A3_nc(2,1) = 0;
A3_nc(2,2) = pi*rho*b^2*U^2;

A3_nc(1,3) = 0;
A3_nc(2,3) = -rho*b^2*U^2*(c*sqrt(1-c^2)-acos(c));
A3_nc(3,1) = 0;
A3_nc(3,2) = -rho*b^2*U^2*(c*sqrt(1-c^2)-acos(c)); % corrected in 2018-01-19
A3_nc(3,3) = -rho*b^2*U^2/pi*(2*c*sqrt(1-c^2)*acos(c)-(1-c^2)-(acos(c))^2);
     
A3_c(1,1) = 0;
A3_c(1,2) = -(2*pi*rho*U*b*U);
A3_c(2,1) = 0;
A3_c(2,2) = -2*pi*rho*U*b^2*1/2*U + b*(a+1/2)*(2*pi*rho*U*b*U);

A3_c(1,3) = -(2*pi*rho*U*b*U/pi*(sqrt(1-c^2)+acos(c)));
A3_c(2,3) = -2*pi*rho*U*b^2*(1/2-(a+1/2))*U/pi*(sqrt(1-c^2)+acos(c));
A3_c(3,1) = 0;
A3_c(3,2) = fat1*U;
A3_c(3,3) = fat1*U/pi*(sqrt(1-c^2)+acos(c));

A4(1,1) = -(2*pi*rho*U*b*1);
A4(1,2) = -(2*pi*rho*U*b*1);
A4(2,1) = b*(a+1/2)*2*pi*rho*U*b*1;
A4(2,2) = b*(a+1/2)*2*pi*rho*U*b*1;
fat2 = -2*rho*U*b^2*(((1+c/2)*sqrt(1-c^2)-(c+1/2)*acos(c))); % common factor 1
A4(3,1) = fat2;
A4(3,2) = fat2;
    
A1 = A1_nc + A1_c;
A2 = A2_nc + A2_c;
A3 = A3_nc + A3_c;

p1 = -0.041*U/b;
q1 = -0.165;
p2 = -0.32*U/b;
q2 = -0.335;

B1(1,1) = q1;
B1(1,2) = q1*b*(1/2-a);
B1(2,1) = q2;
B1(2,2) = q2*b*(1/2-a);

B1(1,3) = q1*b/(2*pi)*((1-2*c)*acos(c)+(2-c)*sqrt(1-c^2));
B1(2,3) = q2*b/(2*pi)*((1-2*c)*acos(c)+(2-c)*sqrt(1-c^2));
  
B2(1,1) = 0;
B2(1,2) = q1*U;
B2(2,1) = 0;
B2(2,2) = q2*U;

B2(1,3) = q1*U/pi*(sqrt(1-c^2)+acos(c));
B2(2,3) = q2*U/pi*(sqrt(1-c^2)+acos(c));

B3(1,1) = 0;
B3(1,2) = 0;
B3(2,1) = 0;
B3(2,2) = 0;

B3(1,3) = 0;
B3(2,3) = 0;
  
B4(1,1) = p1;
B4(1,2) = 0;
B4(2,1) = 0;   
B4(2,2) = p2;




amat.A1 = A1;
amat.A2 = A2;
amat.A3 = A3;
amat.A4 = A4;

amat.A1_nc = A1_nc;
amat.A2_nc = A2_nc;
amat.A3_nc = A3_nc;


amat.B1 = B1;
amat.B2 = B2;
amat.B3 = B3;
amat.B4 = B4;