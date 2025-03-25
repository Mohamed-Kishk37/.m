clc;
clear;
c=299792458;
f0=1;%2.3e9;
la=c/f0;  %wave length

fdev=200e6; %fm chirp diviation;
TCPI=100e-3; %total time of one sweep segment
%to be mesured in the real world appplication
fB1=0.3; %beat frequency 1 in Hz
fB2=-0.2; %beat frequency 2 in Hz

%matrix constants
a=2*fdev/(c*TCPI);
b=2/la;

%coeffecient matrix
A=[a -b; -a -b];

%right hand side (coulomn) vector
B=[fB1;fB2];

%solve linear equation
% x(1)=A(1)/B(1);
% x(2)=A(2)/B(2);

x=A\B;

%unknown vector
%x=[R;vr];

R=x(1);
vr=x(2);

disp(['range [m]: ' num2str(R)]);
disp(['radial velocity [m/sec] : ' num2str(vr)]);
disp(['radial velocity [km/h] : ' num2str(3.6*vr)]);























