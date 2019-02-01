%HW01test
% 22nd Jan 2019. 
% Code to try transforming from joint angles to segment angles. BMEN 619.24
% HW01. 
clear all;clc;

% Jacobian matrix solved on board
J=[-1 0 0 0;-1 1 0 0; 0 1 -1 0; 0 0 -1 1];

Jinv=inv(J);% finding matrix inverse

% Our test values are: 
%Theta=[45; -90; -90; -90];
%Phi=[45; -45 ;45 ;-45];
% We try to go from Theta to Phi, and from Phi to Theta
%% Finding Theta from Phi
Phi=[45; -45 ;45 ;-45];
%Phi=[90;135;-45;5];
Phi2theta=[90;0;0;0];% Offset matrix
Theta_test=J*Phi+Phi2theta;

%% Finding Phi from Theta

Theta=[45; -90; -90; -90];
Theta2phi=[90;90;90;90];% offset matrix
Phi_test=Jinv*Theta+Theta2phi;

% Jinv2=[-1 0 0 0;-1 1 0 0;1 -1 1 0; -1 1 -1 1];
% Theta2phi2=[90; 90;-90; 90];
% Phi_test2=Jinv2*Theta+Theta2phi2;

